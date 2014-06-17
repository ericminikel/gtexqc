#!/bin/bash

gatkjar=/humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar
b37ref=/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta

zcat $wgsvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | sed 's/-[0-9]*$//' | tail -n +10 | sort > wgs.sorted
cat variant.count.outliers | sed 's/-[0-9]*$//' | sort > outliers.sorted
zcat $arrayvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 | sed 's/-[0-9]*-SM.*$//' | sort > array.sorted
cat $exchipvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 | sed 's/-[0-9]*-SM.*$//' | sort > exchip.sorted
cat $wesmeta | tail -n +2 | awk '{print $4}' | sed 's/-[0-9]*$//' | sort | uniq -u > wes.sorted
cat $wgsmeta | tail -n +2 | awk '{print $4}' | sed 's/-[0-9]*$//' | sort | uniq -u > wgs.meta.sorted
cat $wesflag | cut -f1 | sed 's/-[0-9]*-SM.*$//' | sort > flagged.sorted

comm -23 wes.sorted flagged.sorted > keep1
comm -23 keep1 outliers.sorted > keep2
comm -12 keep2 wgs.sorted > keep3
comm -12 keep3 array.sorted > keep4
comm -12 keep4 exchip.sorted > keep5
comm -12 keep5 wgs.meta.sorted > keep6

wc -l keep6 # 140

grep -f keep6 $wesmeta | cut -f38 | uniq -c
    # 140 Agilent
cat $wesmeta | tail -n +2 | awk '{print $4}' | sed 's/-[0-9]*$//' | sort | uniq -d > wes.dups
grep -f wes.dups $wesmeta | cut -f38
# 7 samples which have both ICE and Agilent

grep -f keep6 $wgsmeta | cut -f38 | uniq -c
     # 68 HiSeq 2000
     # 72 HiSeq X

# generate lists of samples to include for different analyses.
comm -12 wes.dups wgs.meta.sorted > alldata.samples # N = 6 samples that have Agilent, ICE, X Ten, 2.5M & ExomeChip
grep -o -f keep6 $wgsmeta | sort | uniq > wgs.samples # N = 68 vs. 72 additional samples with just WGS for 2000 vs. X Ten analysis
cat alldata.samples >> wgs.samples # use the alldata ones for their WGS BAMs too

# BAMs to submit for coverage.
# desired columns:
# 1. Participant ID ($25 = Collaborator ID)
# 2. Longer ID of BAM with everything in it ($1 = Quick Search String (details))
# 3. WES or WGS (hard coded)
# 4. Technology: Agilent, ICE, HiSeq 2000, HiSeq X Ten ($38 = technology)
# 5. BAM path ($31 = BAM path with version)

# grab N = 6*2 exomes, on both Agilent and ICE
grep -f alldata.samples $wesmeta | awk -F"\t" '{print $25"\t"$1"\tWES\t"$38"\t"$31}' > bam.metadata
# grab N = 140 whole genomes, on HiSeq 2000 or X Ten
grep -f wgs.samples $wgsmeta | awk -F"\t" '{print $25"\t"$1"\tWGS\t"$38"\t"$31}' >> bam.metadata
wc -l bam.metadata # 152

cat bam.metadata | cut -f5 | noexist > bams.to.update
grep -v -f bams.to.update bam.metadata > bam.metadata.fixed
grep -f bams.to.update bam.metadata | sed 's/v1/v2/'  >> bam.metadata.fixed

cat bam.metadata.fixed | cut -f5 | noexist # none
cat bam.metadata.fixed | wc -l # 158

# create symlinks with highly unique IDs as BAM names
mkdir -p bams
cat bam.metadata.fixed | awk -F"\t" '{print "ln -s "$5" bams/"$2".bam"}' | bash
ls bams | sed 's/@//' | awk -F"\t" '{print "/humgen/atgu1/fs03/eminikel/053gtexqc/bams/"$0}' > bam.list

# create symlinks to BAI indices of BAMs
cd bams
ll | awk '{print "ln -s "$10" "$8}' | sed 's/bam/bai/g' | bash
cd ..

# get Gencode annotations
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz
cat gencode.v19.annotation.gtf | awk '$3=="CDS" {print $1"\t"$4-1"\t"$5}' | sed 's/^chr//' | sed 's/^M/MT/' > gencode_cds.raw.bed
bedtools sort -i gencode_cds.raw.bed > gencode_cds.sorted.bed
bedtools merge -i gencode_cds.sorted.bed > gencode_cds.bed

# check total length
cat gencode_cds.bed | awk '{sum += $3-$2} END {print sum}' # 35345952
# check number of entries per chrom
cat gencode_cds.bed | cut -f1 | uniq -c

broadexome=/humgen/gsa-hpprojects/GATK/bundle/2.8/b37/Broad.human.exome.b37.interval_list 
cat $broadexome | grep -v ^@ | cut -f1,2,3 > broadexome.raw.bed
bedtools sort -i broadexome.raw.bed > broadexome.sorted.bed
bedtools merge -i broadexome.sorted.bed > broadexome.bed

# check total length
cat broadexome.bed | awk '{sum += $3-$2} END {print sum}' # 32760394
# check number of entries per chrom
cat broadexome.bed | cut -f1 | uniq -c

# make a short interval list for test scripts
cat gencode_cds.bed | head -1000 > test_targets.bed

mkdir -p jobtemp

# just double-check
cat bam.list | noexist
# none! great.


# attempts to run all bams together failed due to lack of memory
# for targetset in {broadexome.bed,gencode_cds.bed}
# do
#     bsub -q priority -P $RANDOM -J gtexqc -M 16000000 \
#         -o jobtemp/job.$targetset.noq.out \
#         -e jobtemp/job.$targetset.noq.err \
#         "java -Xmx16g -jar $gatkjar \
#              -R $b37ref \
#              -T DepthOfCoverage \
#              -o cov_${targetset}_20_1 \
#              -I bam.list \
#              -L $targetset \
#              --omitDepthOutputAtEachBase \
#              --minBaseQuality 20 \
#              --minMappingQuality 1"
#     bsub -q priority -P $RANDOM -J gtexqc -M 16000000 \
#         -o jobtemp/job.$targetset.wq.out \
#         -e jobtemp/job.$targetset.wq.err \
#         "java -Xmx16g -jar $gatkjar \
#              -R $b37ref \
#              -T DepthOfCoverage \
#              -o cov_${targetset}_0_0 \
#              -I bam.list \
#              -L $targetset \
#              --omitDepthOutputAtEachBase \
#              --minBaseQuality 0 \
#              --minMappingQuality 0"
# done

# to be absolutely sure
mkdir -p bybam
for targetset in {broadexome.bed,gencode_cds.bed}
do
    for bam in `cat bam.list | tail -n +2`
    do
        sname=`echo $bam | sed 's/.*\///'`
    bsub -q priority -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/job.bybam.$targetset.$sname.noq.out \
        -e jobtemp/job.bybam.$targetset.$sname.noq.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T DepthOfCoverage \
             -o bybam/cov_${targetset}_${sname}_20_1 \
             -I $bam \
             -L $targetset \
             --omitDepthOutputAtEachBase \
             --minBaseQuality 20 \
             --minMappingQuality 1"
    bsub -q priority -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/job.$targetset.$sname.wq.out \
        -e jobtemp/job.$targetset.$sname.wq.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T DepthOfCoverage \
             -o bybam/cov_${targetset}_${sname}_0_0 \
             -I $bam \
             -L $targetset \
             --omitDepthOutputAtEachBase \
             --minBaseQuality 0 \
             --minMappingQuality 0"
    done
done

# run this block of code every time you need to re-run failed jobs
grep exit jobtemp/*.out | grep 134 | sed 's/\.out:Exited.*//' > redos.txt
mkdir -p jobtemp_failed
for line in `cat redos.txt`
do
    echo $line
    cmd=`grep -A 1 LSBATCH $line.out | tail -n +2`
    bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o $line.out \
        -e $line.err "$cmd"
    mv $line.* jobtemp_failed
done

# the jobs on stuck nodes don't give an informative .out file, need a different method
bjobs -dl > bjobs.dl.20140611.1035.txt
cat bjobs.dl.20140611.1035.txt | sed 's/^ *//' | sed ':a;N;$!ba;s/\n//g' | sed 's/----*/\n/g'| grep Termination | egrep -o "Command <.*>" | sed 's/Command <//' | sed 's/>.*//' > commands-to-resubmit.txt
# needs some further editing due to a few errors
cat commands-to-resubmit.txt | sed 's/Xmx8g/Xmx8g /' | sed 's/omitDepthOutputAtEachBase/omitDepthOutputAtEachBase /' > commands-fixed.txt
while read cmd
do
    outputfile=`echo $cmd | sed 's/.*-o //' | sed 's/-I.*//'`
    outfile=`echo $outputfile".out" | sed 's/bybam/jobtemp/' | sed 's/ //'`
    errfile=`echo $outputfile".err" | sed 's/bybam/jobtemp/' | sed 's/ //'` 
    bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o $outfile \
        -e $errfile \
        "$cmd"
done < commands-fixed.txt

bjobs -q bhour -l > bjobs.l.20140612.0951.txt
cat bjobs.l.20140612.0951.txt  | sed 's/^ *//' | sed ':a;N;$!ba;s/\n//g' | sed 's/----*/\n/g'| grep UNKWN | egrep -o "Command <.*>" | sed 's/Command <//' | sed 's/>.*//' > commands-to-resubmit.txt

zcat $exomesnpvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 > exome.snp.cols
zcat $wgsvcf      | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 > wgs.cols
zcat $exomesuppsnpvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 > exome.supp.snp.cols
zcat $arrayvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 > array.cols
cat $exchipvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 > exchip.cols

grep -f alldata.samples exome.snp.cols > exome.ice.ids
grep -f alldata.samples exome.supp.snp.cols > exome.agilent.ids
grep -f wgs.samples wgs.cols > wgs.ids
grep -f wgs.samples array.cols > array.ids
grep -f wgs.samples exchip.cols > exchip.ids

cat $wgsmeta | grep "HiSeq 2000" | grep -o -f wgs.ids - | sort | uniq > wgs.2000.ids
cat $wgsmeta | grep "HiSeq X" | grep -o -f wgs.ids - | sort | uniq > wgs.xten.ids


# WGS
bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/wgs.2000.snp.out \
        -e jobtemp/wgs.2000.snp.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType SNP \
             -sf wgs.2000.ids \
             -env \
             -o wgs.2000.snp.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/wgs.xten.snp.out \
        -e jobtemp/wgs.xten.snp.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType SNP \
             -sf wgs.xten.ids \
             -env \
             -o wgs.xten.snp.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/wgs.2000.indel.out \
        -e jobtemp/wgs.2000.indel.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType INDEL \
             -sf wgs.2000.ids \
             -env \
             -o wgs.2000.indel.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/wgs.xten.indel.out \
        -e jobtemp/wgs.xten.indel.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType INDEL \
             -sf wgs.xten.ids \
             -env \
             -o wgs.xten.indel.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/exome.ice.snp.out \
        -e jobtemp/exome.ice.snp.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomesnpvcf \
             -sf exome.ice.ids \
             -env \
             -o exome.ice.snp.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/exome.ice.indel.out \
        -e jobtemp/exome.ice.indel.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomeindelvcf \
             -sf exome.ice.ids \
             -env \
             -o exome.ice.indel.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/exome.agilent.snp.out \
        -e jobtemp/exome.agilent.snp.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomesuppsnpvcf \
             -sf exome.agilent.ids \
             -env \
             -o exome.agilent.snp.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/exome.agilent.indel.out \
        -e jobtemp/exome.agilent.indel.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomesuppindelvcf \
             -sf exome.agilent.ids \
             -env \
             -o exome.agilent.indel.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/array.snp.out \
        -e jobtemp/array.snp.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $arrayvcf \
             -sf array.ids \
             -env \
             -o array.vcf"

# bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
#         -o jobtemp/exchip.snp.out \
#         -e jobtemp/exchip.snp.err \
# "java -Xmx8g -jar $gatkjar \
#              -R $b37ref \
#              -T SelectVariants \
#              -V $exchipvcf \
#              -sf exchip.ids \
#              -o exchip.vcf"

##### ERROR MESSAGE: Key GCR found in VariantContext field INFO at 1:564766 but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.
cat $exchipvcf | sed 's/GCR=[0-9\.]*;//' > exchip.no.gcr.vcf # takes only a few seconds
##### ERROR MESSAGE: Key MN found in VariantContext field FILTER at 1:569406 but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.
sed -i 's/\tMN\t/\t.\t/' exchip.no.gcr.vcf
##### ERROR MESSAGE: Key hh found in VariantContext field INFO at X:1404785 but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.
cat $exchipvcf | grep ^# > exchip.fixed.vcf # preserve header
# now replace FILTER and INFO with . and nothing
cat $exchipvcf | grep -v ^# | awk -F"\t" 'BEGIN {OFS=FS} $7="."; ' | awk -F"\t" 'BEGIN {OFS=FS} $8=".";' >> exchip.fixed.vcf
bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
         -o jobtemp/exchip.fixed.snp.out \
         -e jobtemp/exchip.fixed.snp.err \
 "java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T SelectVariants \
              -V exchip.fixed.vcf \
              -sf exchip.ids \
              -env \
              -o exchip.vcf"

# monkol's suggestion:
# -U LENIENT_VCF_PROCESSING \
# -S LENIENT \




### REPLACING OF SAMPLE IDs

# figure out how to replace sample ids
# cat $exchipvcf | grep -m 1 ^#CHROM | sed 's/-[0-9]*-SM-[0-9A-Z]*//g'
# cat $exchipvcf | grep ^# > exchipvcfheader.vcf
# sed -i 's/-[0-9]*-SM-[0-9A-Z]*//g' exchipvcfheader.vcf # does the replacement in every line
# cat $exchipvcf | grep ^# > exchipvcfheader.vcf
# # sed -i '2,3s/-[0-9]*-SM-[0-9A-Z]*//g' exchipvcfheader.vcf
# chromlineno=`grep -n ^#CHROM exchipvcfheader.vcf | cut -f1 -d:`
# sed ${chromlineno}'s/-[0-9]*-SM-[0-9A-Z]*//g' exchipvcfheader.vcf

# ok, now actually replace all the sample ids
# array: remove strings like -0002-SM-3PYKN 
chromlineno=`grep -m 1 -n ^#CHROM array.vcf | cut -f1 -d:`
head -n $chromlineno array.vcf | sed ${chromlineno}'s/-[0-9]*-SM-[0-9A-Z]*//g' > header.txt
bgzip -c array.vcf > array.vcf.bgz
tabix -r header.txt array.vcf.bgz > array.sn.vcf.bgz

# exchip: remove strings like -0002-SM-3PYKN 
chromlineno=`grep -m 1 -n ^#CHROM exchip.vcf | cut -f1 -d:`
head -n $chromlineno exchip.vcf | sed ${chromlineno}'s/-[0-9]*-SM-[0-9A-Z]*//g' > header.txt
bgzip -c exchip.vcf > exchip.vcf.bgz
tabix -r header.txt exchip.vcf.bgz > exchip.sn.vcf.bgz
#sed ${chromlineno}'s/-[0-9]*-SM-[0-9A-Z]*//g' exchip.vcf > exchip.sn.vcf


# gunzip exome.*.vcf.gz
for exomefile in exome.*.vcf
do
    # exome: remove C1422:: or C1587:: and strings like -0002
    chromlineno=`grep -m 1 -n ^#CHROM $exomefile | cut -f1 -d:`
    head -n $chromlineno $exomefile | sed -r ${chromlineno}'s/-[0-9]+//g' | sed ${chromlineno}'s/C[0-9]*:://g' > header.txt
    bgzip $exomefile
    newfilename=`echo $exomefile | sed 's/vcf/sn.vcf.bgz/'`
    tabix -r header.txt $exomefile.gz > $newfilename
done

# wgsfile=wgs.xten.indel.vcf
# for WGS, first bgzip them
for wgsfile in wgs.*.vcf
do
    bsub -q bweek -P $RANDOM -J bgzwgs -M 8000000 \
            -o jobtemp/bgz.$wgsfile.out \
            -e jobtemp/bgz.$wgsfile.err \
    "bgzip $wgsfile"
done
# if you need to re-run and force overwrite, change bgzip to bgzip -f
# then come back and reheader them after all have finished bgzipping
for wgsfile in wgs.*.vcf.gz
do
    # genome: just remove strings like -0002
    bsub -q bweek -P $RANDOM -J snwgs -M 8000000 \
            -o jobtemp/sn.$wgsfile.out \
            -e jobtemp/sn.$wgsfile.err \
    "chromlineno=`zcat $wgsfile | grep -m 1 -n ^#CHROM - | cut -f1 -d:`
    zcat $wgsfile | head -n \$chromlineno | sed -r \${chromlineno}'s/-[0-9]+//g' > header.$wgsfile.txt
    newfilename=`echo $wgsfile | sed 's/vcf.gz/sn.vcf.bgz/'`
    tabix -r header.$wgsfile.txt $wgsfile > \$newfilename"
done
# these jobs took 87 to 623 seconds.

# this one was used as a test
# wgsfile=scratch.vcf.gz
#     bsub -q bweek -P $RANDOM -J snwgs -M 8000000 \
#             -o /dev/null \
#             -e /dev/null \
#     "chromlineno=`zcat $wgsfile | grep -m 1 -n ^#CHROM - | cut -f1 -d:`
#     zcat $wgsfile | head -n \$chromlineno | sed -r \${chromlineno}'s/-[0-9]+//g' > header.$wgsfile.txt
#     newfilename=`echo $wgsfile | sed 's/vcf.gz/sn.vcf.bgz/'`
#     tabix -r header.$wgsfile.txt $wgsfile > \$newfilename"

# re-casting of the bybam coverage data
# only need 3 things:
# - coverage_proportions
# - coverage_counts
# - interval_summary - only the mean coverage stats
cd bybam
# coverage proportions - one row per sample
fname=`ls *coverage_proportions | head -1`
echo -n -e "sid\ttargets\tqual\t" > coverage_proportions.txt
cat $fname | head -1 | cut -f 2- >> coverage_proportions.txt
for fname in *coverage_proportions
do
    sid=`echo $fname | sed 's/cov_[a-z_]*.bed_//' | sed 's/\.bam.*//'` # long sample id
    targets=`echo $fname | sed 's/cov_//' | sed 's/\.bed.*//'` # gencode_cds or broadexome
    qual=`echo $fname | sed 's/.*bam_//' | sed 's/\.sample.*//'` # 20_1 or 0_0
    echo -n -e "$sid\t$targets\t$qual\t" >> coverage_proportions.txt
    cat $fname | tail -1 | cut -f 2- >> coverage_proportions.txt
done

# coverage counts - one row per sample; exactly like coverage_proportions above.
fname=`ls *coverage_counts | head -1`
echo -n -e "sid\ttargets\tqual\t" > coverage_counts.txt
cat $fname | head -1 | cut -f 2- >> coverage_counts.txt
for fname in *coverage_counts
do
    sid=`echo $fname | sed 's/cov_[a-z_]*.bed_//' | sed 's/\.bam.*//'` # long sample id
    targets=`echo $fname | sed 's/cov_//' | sed 's/\.bed.*//'` # gencode_cds or broadexome
    qual=`echo $fname | sed 's/.*bam_//' | sed 's/\.sample.*//'` # 20_1 or 0_0
    echo -n -e "$sid\t$targets\t$qual\t" >> coverage_counts.txt
    cat $fname | tail -1 | cut -f 2- >> coverage_counts.txt
done

# for interval_summary, we'll want to transpose, and handle each targetset separately
for targetset in {broadexome.bed,gencode_cds.bed}
do
    fname=`ls *$targetset*interval_summary | head -1`
    echo -n -e "sid\tqual\t" > interval_summary_$targetset.txt
    cat $fname | cut -f1 | tail -n +2 | tr '\n' '\t' >> interval_summary_$targetset.txt
    echo -n -e "\n" >> interval_summary_$targetset.txt
    for fname in *$targetset*interval_summary
    do
        sid=`echo $fname | sed 's/cov_[a-z_]*.bed_//' | sed 's/\.bam.*//'` # long sample id
        qual=`echo $fname | sed 's/.*bam_//' | sed 's/\.sample.*//'` # 20_1 or 0_0
        echo -n -e "$sid\t$qual\t" >> interval_summary_$targetset.txt 
        cat $fname | cut -f5 | tail -n +2 | tr '\n' '\t' >> interval_summary_$targetset.txt # $5 is sample avg coverage on interval
        echo -n -e "\n" >> interval_summary_$targetset.txt
    done
done
cd ..

# demonstrate that .bgz files must be renamed to .gz for GATK to accept them
java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T ValidateVariants \
   -V exchip.sn.vcf.bgz

mv exchip.sn.vcf.bgz exchip.sn.vcf.gz
tabix exchip.sn.vcf.gz

java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T ValidateVariants \
   -V exchip.sn.vcf.gz

# now do this for all bgz files
for fname in *.bgz
do
    newfname=`echo $fname | sed 's/bgz/gz/'`
    mv $fname $newfname
done

for fname in *sn.vcf.gz
do
    bsub -q bhour -P $RANDOM -J tabix -M 8000000 \
            -o jobtemp/tabix.$fname.out \
            -e jobtemp/tabix.$fname.err \
    "tabix -f $fname"
done

#### genotypic concordance

# compare each sequencing dataset with the 2.5M array as "truth"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.wgs.xten.vs.array.gq30dp10.molt.out \
            -e jobtemp/conc.wgs.xten.vs.array.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp array.sn.vcf.gz \
              -eval wgs.xten.snp.sn.vcf.gz \
              -moltenize \
              -o wgs.xten.vs.array.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.wgs.2000.vs.array.gq30dp10.molt.out \
            -e jobtemp/conc.wgs.2000.vs.array.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp array.sn.vcf.gz \
              -eval wgs.2000.snp.sn.vcf.gz \
              -moltenize \
              -o wgs.2000.vs.array.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.exome.agilent.vs.array.gq30dp10.molt.out \
            -e jobtemp/conc.exome.agilent.vs.array.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp array.sn.vcf.gz \
              -eval exome.agilent.snp.sn.vcf.gz \
              -moltenize \
              -o exome.agilent.vs.array.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.exome.ice.vs.array.gq30dp10.molt.out \
            -e jobtemp/conc.exome.ice.vs.array.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp array.sn.vcf.gz \
              -eval exome.ice.snp.sn.vcf.gz \
              -moltenize \
              -o exome.ice.vs.array.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.exchip.vs.array.gq30dp10.molt.out \
            -e jobtemp/conc.exchip.vs.array.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp array.sn.vcf.gz \
              -eval exchip.sn.vcf.gz \
              -moltenize \
              -o exchip.vs.array.gq30dp10.molt"

# create an interval list of only the 2.5M array SNP sites
zcat array.sn.vcf.gz | grep -v ^# | head | cut -f1-5 # work in progress


# figure out which samples dropped from the 2.5M array on QC
cat $arrayfam_postqc | sort | cut -d ' ' -f1 > postqc.indivs
cat $arrayfam_preqc  | sort | cut -d ' ' -f1 > preqc.indivs
comm -23 preqc.indivs postqc.indivs # gives 4 differences

# are any of these in the current analysis set?
comm -23 preqc.indivs postqc.indivs | sed 's/-[0-9]*-SM-.*//' > array.qc.removed.indivs
grep -f array.qc.removed.indivs alldata.samples # none, good
grep -f array.qc.removed.indivs wgs.samples # none, good


# entire interval_summary table appears too large to read into R
# therefore attempt to separate out the interval_summary data by technology
# and then take averages using command line tools
# code to take average of each column in awk from: http://www.unix.com/shell-programming-and-scripting/186493-awk-based-script-find-average-all-columns-data-file.html

for targetset in {broadexome.bed,gencode_cds.bed}
do
    # separate out the high and low quality depth calculations
    cat bybam/interval_summary_${targetset}.txt | awk '$2 == "20_1" {print}' > bybam/interval_summary_${targetset}_hq.txt
    cat bybam/interval_summary_${targetset}.txt | awk '$2 == "0_0" {print}' > bybam/interval_summary_${targetset}_lq.txt
    
    # create greppable (^ at front of line) files listing samples of interest by technology
    cat $wesmeta | grep ICE | cut -f1 | grep -f alldata.samples - | awk '{print "^"$1}' > ice.sids.grepready
    cat $wesmeta | grep Agilent | cut -f1 | grep -f alldata.samples - | awk '{print "^"$1}' > agilent.sids.grepready
    cat $wgsmeta | grep "HiSeq 2000" | cut -f1 | grep -f wgs.samples - | awk '{print "^"$1}' > h2.sids.grepready
    cat $wgsmeta | grep "HiSeq X" | cut -f1 | grep -f wgs.samples - | awk '{print "^"$1}' > hx.sids.grepready
    
    # extract sample sets by quality and technology
    grep -f ice.sids.grepready bybam/interval_summary_${targetset}_hq.txt > bybam/interval_summary_${targetset}_hq_ice.txt
    grep -f ice.sids.grepready bybam/interval_summary_${targetset}_lq.txt > bybam/interval_summary_${targetset}_lq_ice.txt
    grep -f agilent.sids.grepready bybam/interval_summary_${targetset}_hq.txt > bybam/interval_summary_${targetset}_hq_agilent.txt
    grep -f agilent.sids.grepready bybam/interval_summary_${targetset}_lq.txt > bybam/interval_summary_${targetset}_lq_agilent.txt
    grep -f h2.sids.grepready bybam/interval_summary_${targetset}_hq.txt > bybam/interval_summary_${targetset}_hq_h2.txt
    grep -f h2.sids.grepready bybam/interval_summary_${targetset}_lq.txt > bybam/interval_summary_${targetset}_lq_h2.txt
    grep -f hx.sids.grepready bybam/interval_summary_${targetset}_hq.txt > bybam/interval_summary_${targetset}_hq_hx.txt
    grep -f hx.sids.grepready bybam/interval_summary_${targetset}_lq.txt > bybam/interval_summary_${targetset}_lq_hx.txt
    
    # take column means
    cat bybam/interval_summary_${targetset}_hq_ice.txt     | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_hq_ic_means.txt
    cat bybam/interval_summary_${targetset}_lq_ice.txt     | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_lq_ic_means.txt
    cat bybam/interval_summary_${targetset}_hq_agilent.txt | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_hq_ag_means.txt
    cat bybam/interval_summary_${targetset}_lq_agilent.txt | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_lq_ag_means.txt
    cat bybam/interval_summary_${targetset}_hq_h2.txt      | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_hq_h2_means.txt
    cat bybam/interval_summary_${targetset}_lq_h2.txt      | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_lq_h2_means.txt
    cat bybam/interval_summary_${targetset}_hq_hx.txt      | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_hq_hx_means.txt
    cat bybam/interval_summary_${targetset}_lq_hx.txt      | awk '{f=NF;for(i=3;i<=NF;i++)a[i]+=$i}END{for(i=3;i<=f;i++)printf a[i]/(NR)"\t";print ""}' > is_${targetset}_lq_hx_means.txt
    
    # transpose
    cat is_${targetset}_hq_ic_means.txt | tr '\t' '\n' > is_${targetset}_hq_ic_means_t.txt
    cat is_${targetset}_lq_ic_means.txt | tr '\t' '\n' > is_${targetset}_lq_ic_means_t.txt
    cat is_${targetset}_hq_ag_means.txt | tr '\t' '\n' > is_${targetset}_hq_ag_means_t.txt
    cat is_${targetset}_lq_ag_means.txt | tr '\t' '\n' > is_${targetset}_lq_ag_means_t.txt
    cat is_${targetset}_hq_h2_means.txt | tr '\t' '\n' > is_${targetset}_hq_h2_means_t.txt
    cat is_${targetset}_lq_h2_means.txt | tr '\t' '\n' > is_${targetset}_lq_h2_means_t.txt
    cat is_${targetset}_hq_hx_means.txt | tr '\t' '\n' > is_${targetset}_hq_hx_means_t.txt
    cat is_${targetset}_lq_hx_means.txt | tr '\t' '\n' > is_${targetset}_lq_hx_means_t.txt
    
    # check
    wc -l is_${targetset}_hq_ic_means_t.txt
    wc -l is_${targetset}_lq_ic_means_t.txt
    wc -l is_${targetset}_hq_ag_means_t.txt
    wc -l is_${targetset}_lq_ag_means_t.txt
    wc -l is_${targetset}_hq_h2_means_t.txt
    wc -l is_${targetset}_lq_h2_means_t.txt
    wc -l is_${targetset}_hq_hx_means_t.txt
    wc -l is_${targetset}_lq_hx_means_t.txt
    # ok
    
    # paste together, with an appropriate header and first column
    header_source=`ls bybam/*${targetset}*interval_summary | head -1`
    cat $header_source | cut -f1 | tail -n +2 > is_${targetset}_interval_names.txt
    echo -ne "\n" >> is_${targetset}_interval_names.txt # add one blank row to match the data
    wc -l is_${targetset}_interval_names.txt
    echo -e "interval\tICE.20_1\tICE.0_0\tAgilent.20_1\tAgilent.0_0\th2000.20_1\th2000.0_0\thX.20_1\thX.0_0" > is_${targetset}_all_means.txt
    paste is_${targetset}_interval_names.txt is_${targetset}_hq_ic_means_t.txt is_${targetset}_lq_ic_means_t.txt is_${targetset}_hq_ag_means_t.txt is_${targetset}_lq_ag_means_t.txt is_${targetset}_hq_h2_means_t.txt is_${targetset}_lq_h2_means_t.txt is_${targetset}_hq_hx_means_t.txt is_${targetset}_lq_hx_means_t.txt >> is_${targetset}_all_means.txt
done

#### check on DepthOfCoverage jobs
# 2014-06-11 09:65 check-in
# 158*4 = 632 files should be created
ls bybam/*coverage_counts | wc -l #606
bjobs | grep -v JOBID | grep -v interact | wc -l #26
# 606 + 26 = 632, all accounted for.

ls bybam/*coverage_counts | wc -l
# 632 # huzzah - all finished as of June 12, 1:37p

############

#### get summary stats from concordance:
for fname in *.molt
do
    cat $fname | grep -v GTEX > $fname.summary
done

fname=`ls *.molt.summary | head -1`
echo -en "comparison\t" > all.concordance.summary
cat $fname | grep -A 2 ^#:GATKTable:4 | tail -1 >> all.concordance.summary
for fname in *.molt.summary
do
    echo -en $fname"\t" >> all.concordance.summary
    cat $fname | grep -A 3 ^#:GATKTable:4 | tail -1 >> all.concordance.summary
done

for fname in *.molt.summary
do
    cat $fname | head -23 > $fname.concordance.proportions
done

# do genotype concordance against 2.5M array with target sites only
for seqset in {wgs.xten,wgs.2000,exome.agilent,exome.ice}
do
    for targetset in {broadexome.bed,gencode_cds.bed}
    do
    bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.${targetset}.${seqset}.vs.array.gq30dp10.molt.out \
            -e jobtemp/conc.${targetset}.${seqset}.vs.array.gq30dp10.molt.err \
            "java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -L $targetset \
              -comp array.sn.vcf.gz \
              -eval ${seqset}.snp.sn.vcf.gz \
              -moltenize \
              -o gc.${seqset}.vs.array.${targetset}.bed.gq30dp10.moltgc"
    done
done

fname=`ls *.moltgc | head -1`
echo -en "comparison\t" > interval.concordance.summary
cat $fname | grep -A 2 ^#:GATKTable:4 | tail -1 >> interval.concordance.summary
for fname in *.moltgc
do
    echo -en $fname"\t" >> interval.concordance.summary
    cat $fname | grep -A 3 ^#:GATKTable:4 | tail -1 >> interval.concordance.summary
done



# check sex of people in diff. data subsets
# use sex in $arrayfam_postqc
cat bam.metadata.fixed | grep "HiSeq 2000" | cut -f1 > h2000.snames
grep -f h2000.snames $arrayfam_postqc | awk '$5=="1" {print "M"} $5=="2" {print "F"}' | sort | uniq -c
     # 28 F
     # 40 M
cat bam.metadata.fixed | grep "HiSeq X" | cut -f1 > hX.snames
grep -f hX.snames $arrayfam_postqc | awk '$5=="1" {print "M"} $5=="2" {print "F"}' | sort | uniq -c
     # 22 F
     # 56 M
grep -f alldata.samples $arrayfam_postqc | awk '$5=="1" {print "M"} $5=="2" {print "F"}' | sort | uniq -c
      # 1 F
      # 5 M


# final genotypic concordance comparison: ICE vs. X Ten, in the Gencode CDS.
bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.ice.vs.xten.snp.gq30dp10.molt.out \
            -e jobtemp/conc.ice.vs.xten.snp.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp wgs.xten.snp.sn.vcf.gz \
              -eval exome.ice.snp.sn.vcf.gz  \
              -moltenize \
              -o ice.vs.xten.snp.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.ice.vs.xten.indel.gq30dp10.molt.out \
            -e jobtemp/conc.ice.vs.xten.indel.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp wgs.xten.indel.sn.vcf.gz \
              -eval exome.ice.indel.sn.vcf.gz  \
              -moltenize \
              -o ice.vs.xten.indel.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.agilent.vs.xten.snp.gq30dp10.molt.out \
            -e jobtemp/conc.agilent.vs.xten.snp.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp wgs.xten.snp.sn.vcf.gz \
              -eval exome.agilent.snp.sn.vcf.gz  \
              -moltenize \
              -o agilent.vs.xten.snp.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J conc -M 8000000 \
            -o jobtemp/conc.agilent.vs.xten.indel.gq30dp10.molt.out \
            -e jobtemp/conc.agilent.vs.xten.indel.gq30dp10.molt.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp wgs.xten.indel.sn.vcf.gz \
              -eval exome.agilent.indel.sn.vcf.gz  \
              -moltenize \
              -o agilent.vs.xten.indel.gq30dp10.molt"

cat ice.vs.xten.indel.gq30dp10.molt     | grep -A 2 ^#:GATKTable:GenotypeConcordance_Summary
cat agilent.vs.xten.indel.gq30dp10.molt | grep -A 2 ^#:GATKTable:GenotypeConcordance_Summary
cat ice.vs.xten.snp.gq30dp10.molt     | grep -A 2 ^#:GATKTable:GenotypeConcordance_Summary
cat agilent.vs.xten.snp.gq30dp10.molt | grep -A 2 ^#:GATKTable:GenotypeConcordance_Summary



java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -comp array.sn.vcf.gz \
              -eval wgs.xten.snp.sn.vcf.gz \
              -moltenize \
              -o gc.wgs.xten.vs.array.gencode_cds.bed.molt2

cat gc.wgs.xten.vs.array.gencode_cds.bed.molt2 | less

cat gc.wgs.xten.vs.array.gencode_cds.bed.molt2 | head -23 > gc.wgs.xten.vs.array.gencode_cds.bed.molt2.concordance.proportions



#### adding ICE v. XTen analysis with N>6

# 1. find individuals with both ICE and X ten

cat $wesmeta | grep ICE | wc -l
# 351
cat $wgsmeta | grep "HiSeq X" | wc -l
# 85

cat $wesmeta | grep ICE | awk -F"\t" '{print $25}' | sort > ice.exome.snames
cat $wgsmeta | grep "HiSeq X" | awk -F"\t" '{print $25}' | sort > hx.genome.snames
comm -12 ice.exome.snames hx.genome.snames | wc -l
# 13

cat $wesmeta | grep Agilent | awk -F"\t" '{print $25}' | sort > agilent.exome.snames
comm -12 agilent.exome.snames hx.genome.snames > ag-hx.snames
wc -l ag-hx.snames
# 78

# there are 85 HiSeq X Ten genomes. 78 have Agilent, 13 have ICE, 6 have both

# grep -f ag-hx.snames $wesmeta | grep Agilent   | awk -F"\t" '{print $25"\t"$1"\tWES\t"$38"\t"$31}' >  ag-hx.metadata
# grep -f ag-hx.snames $wgsmeta | grep "HiSeq X" | awk -F"\t" '{print $25"\t"$1"\tWES\t"$38"\t"$31}' >> ag-hx.metadata

grep -f ag-hx.snames exome.supp.cols

grep -f alldata.samples $wesmeta | awk -F"\t" '{print $25"\t"$1"\tWES\t"$38"\t"$31}' > bam.metadata
# grab N = 140 whole genomes, on HiSeq 2000 or X Ten
grep -f wgs.samples $wgsmeta | awk -F"\t" '{print $25"\t"$1"\tWGS\t"$38"\t"$31}' >> bam.metadata
wc -l bam.metadata # 152

grep -f ag-hx.snames exome.snp.cols | grep -v 1587 | sort > ag-hx.agilent.cols
grep -f ag-hx.snames wgs.cols | sort > ag-hx.xten.cols
comm -12 ag-hx.xten.cols ag-hx.agilent.cols > ag-hx.use.cols
# 72! awesome, they *already* have the same IDs.

# subset to only these 72 individuals
bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.exome.agilent.snp.out \
        -e jobtemp/ag-hx.exome.agilent.snp.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomesnpvcf \
             -sf ag-hx.use.cols \
             -env \
             -o ag-hx.exome.agilent.snp.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.exome.agilent.indel.out \
        -e jobtemp/ag-hx.exome.agilent.indel.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomeindelvcf \
             -sf ag-hx.use.cols \
             -env \
             -o ag-hx.exome.agilent.indel.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.wgs.xten.indel.out \
        -e jobtemp/ag-hx.wgs.xten.indel.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType INDEL \
             -sf ag-hx.use.cols \
             -env \
             -o ag-hx.wgs.xten.indel.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.wgs.xten.snp.out \
        -e jobtemp/ag-hx.wgs.xten.snp.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType SNP \
             -sf ag-hx.use.cols \
             -env \
             -o ag-hx.wgs.xten.snp.vcf"


# bgzip and tabix the VCFs.
for fname in ag-hx.*snp.vcf
do
    bsub -q bhour -W 4:00 -P $RANDOM -J bgztb -M 8000000 \
        -o jobtemp/$fname.bgzip.tabix.out \
        -e jobtemp/$fname.bgzip.tabix.err \
        "bgzip $fname; tabix $fname.gz"
done




bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.xten.vs.agilent.indel.out \
        -e jobtemp/ag-hx.xten.vs.agilent.indel.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp ag-hx.exome.agilent.indel.vcf.gz \
              -eval ag-hx.wgs.xten.indel.vcf.gz  \
              -moltenize \
              -o ag-hx.xten.vs.agilent.indel.gq30dp10.molt"

# next need tosubmit this:
bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.xten.vs.agilent.snp.out \
        -e jobtemp/ag-hx.xten.vs.agilent.snp.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp ag-hx.exome.agilent.snp.vcf.gz \
              -eval ag-hx.wgs.xten.snp.vcf.gz  \
              -moltenize \
              -o ag-hx.xten.vs.agilent.snp.gq30dp10.molt"


# now switch with X Ten as standard and Agilent being evaluted:
bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.agilent.vs.xten.indel.out \
        -e jobtemp/ag-hx.agilent.vs.xten.indel.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp ag-hx.wgs.xten.indel.vcf.gz \
              -eval ag-hx.exome.agilent.indel.vcf.gz \
              -moltenize \
              -o ag-hx.agilent.vs.xten.indel.gq30dp10.molt"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-hx.agilent.vs.xten.snp.out \
        -e jobtemp/ag-hx.agilent.vs.xten.snp.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp ag-hx.wgs.xten.snp.vcf.gz \
              -eval ag-hx.exome.agilent.snp.vcf.gz \
              -moltenize \
              -o ag-hx.agilent.vs.xten.snp.gq30dp10.molt"


# now same with HiSeq 2000


cat $wgsmeta | grep "HiSeq 2000" | wc -l
# 80
cat $wgsmeta | grep "HiSeq 2000" | awk -F"\t" '{print $25}' | sort > h2.genome.snames
comm -12 agilent.exome.snames h2.genome.snames > ag-h2.snames
# 71

grep -f ag-h2.snames exome.snp.cols | grep -v 1587 | sort > ag-h2.agilent.cols
grep -f ag-h2.snames wgs.cols | sort > ag-h2.2000.cols
comm -12 ag-h2.2000.cols ag-h2.agilent.cols > ag-h2.use.cols1
wc -l ag-h2.use.cols1
# 71. for simplicity let's keep only the 68 from the main analysis
grep -f wgs.samples ag-h2.use.cols1 > ag-h2.use.cols
wc -l ag-h2.use.cols
# 68

# subset to only these 72 individuals
bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-h2.exome.agilent.snp.out \
        -e jobtemp/ag-h2.exome.agilent.snp.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomesnpvcf \
             -sf ag-h2.use.cols \
             -env \
             -o ag-h2.exome.agilent.snp.vcf"

bsub -q bhour -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-h2.exome.agilent.indel.out \
        -e jobtemp/ag-h2.exome.agilent.indel.err \
"java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $exomeindelvcf \
             -sf ag-h2.use.cols \
             -env \
             -o ag-h2.exome.agilent.indel.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-h2.wgs.2000.indel.out \
        -e jobtemp/ag-h2.wgs.2000.indel.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType INDEL \
             -sf ag-h2.use.cols \
             -env \
             -o ag-h2.wgs.2000.indel.vcf"

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-h2.wgs.2000.snp.out \
        -e jobtemp/ag-h2.wgs.2000.snp.err \
        "java -Xmx8g -jar $gatkjar \
             -R $b37ref \
             -T SelectVariants \
             -V $wgsvcf \
             -selectType SNP \
             -sf ag-h2.use.cols \
             -env \
             -o ag-h2.wgs.2000.snp.vcf"

# NEXT DO:

# bgzip and tabix the VCFs.
for fname in ag-h2.*snp.vcf
do
    bsub -q bhour -W 4:00 -P $RANDOM -J bgztb -M 8000000 \
        -o jobtemp/$fname.bgzip.tabix.out \
        -e jobtemp/$fname.bgzip.tabix.err \
        "bgzip $fname; tabix $fname.gz"
done

# AND THEN:

bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-h2.2000.vs.agilent.indel.out \
        -e jobtemp/ag-h2.2000.vs.agilent.indel.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp ag-h2.exome.agilent.indel.vcf.gz \
              -eval ag-h2.wgs.2000.indel.vcf.gz  \
              -moltenize \
              -o ag-h2.2000.vs.agilent.indel.gq30dp10.molt"

# next need tosubmit this:
bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/ag-h2.2000.vs.agilent.snp.out \
        -e jobtemp/ag-h2.2000.vs.agilent.snp.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -L gencode_cds.bed \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp ag-h2.exome.agilent.snp.vcf.gz \
              -eval ag-h2.wgs.2000.snp.vcf.gz  \
              -moltenize \
              -o ag-h2.2000.vs.agilent.snp.gq30dp10.molt"


# GC content of Gencode intervals, to see if this is correlated with
# loss of coverage in X Ten vs. 2000
java -Xmx8g -jar $gatkjar \
   -T GCContentByInterval \
   -R $b37ref \
   -L gencode_cds.bed \
   -o gencode_cds.gccontent.txt

# compare the MIXED sites to the 2.5M array sites
bsub -q bweek -P $RANDOM -J gtexqc -M 8000000 \
        -o jobtemp/wgs.mixed.vcf.out \
        -e jobtemp/wgs.mixed.vcf.err \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T SelectVariants \
              -V $wgsvcf \
              -selectType MIXED \
              -o wgs.mixed.only.vcf"


