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
grep exit jobtemp/*.out | sed 's/\.out:Exited.*//' > redos.txt
mkdir -p jobtemp_failed
for line in `cat redos.txt`
do
    echo $line
    cmd=`grep -A 1 LSBATCH $line.out | tail -n +2`
    bsub -q bhour -W 4:00 -P $RANDOM -J gtexqc -M 8000000 \
        -o $line.out \
        -e $line.err "$cmd"
    mv $line.* jobtemp_failed
done

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

# for WGS, first bgzip them
for wgsfile in wgs.*.vcf
do
    bsub -q bweek -P $RANDOM -J bgzwgs -M 8000000 \
            -o jobtemp/bgz.$wgsfile.out \
            -e jobtemp/bgz.$wgsfile.err \
    "bgzip $wgsfile"
done
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




