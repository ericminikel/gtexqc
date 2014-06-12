
#### prefatory remarks
setwd("~/d/sci/053gtexqc/data")
options(stringsAsFactors=FALSE)
require(sqldf)
require(reshape2)
require(ggplot2)

midpoints = function(numericvec) {
  midpointvec = numeric(length(numericvec)-1)
  midpointvec = 0.5*(numericvec[1:(length(numericvec)-1)]+
                     numericvec[2:length(numericvec)])
  return(midpointvec)
}


#### style parameters
gcolor = '#FFA824' # aureoline yellow
ecolor = '#0D4F8B' # indigo dye

newlty = 1 # solid line for the newer of the technologies
oldlty = 3 # dotted line for the older of the technologies

hqcolor = '#000000' # black for high-quality coverage
lqcolor = '#999999' # gray for all coverage incl low quality

getk = function(ge) {
  if (ge=="WES") {
    return (ecolor)
  } else if (ge=="WGS") {
    return (gcolor)
  } else {
    stop("invalid specification: should be WES or WGS")
  }
}

getlty = function(tech) {
  if (tech %in% c("ICE","HiSeq X")) {
    return (newlty)
  } else if (tech %in% c("Agilent","HiSeq 2000")) {
    return (oldlty)
  } else {
    stop("invalid specification: should be Agilent, ICE, HiSeq 2000 or HiSeq X")
  }
}

legtable = data.frame(
  text=c("ICE Exome","Agilent Exome","HiSeq X Ten Genome","HiSeq 2000 Genome"),
  lty=c(newlty,oldlty,newlty,oldlty),
  k=c(ecolor,ecolor,gcolor,gcolor))
# this table allows you to do things like this:
# legend("bottomleft",legtable$text,col=legtable$k,lty=legtable$lty,lwd=2)

# named vector of full names for plot titles and stuff. like a Python dictionary.
fullname = c("Gencode coding intervals","Broad Exome intervals","whole exomes","whole genomes",
             "ICE exomes","Agilent exomes","HiSeq 2000 genomes","HiSeq X Ten genomes",
             "BQ ≥ 20 & MAPQ ≥ 1","No quality filters")
names(fullname) = c('gencode_cds','broadexome','WES','WGS','ICE','Agilent','HiSeq 2000','HiSeq X',
                    '20_1','0_0')


#### read in data
cp = read.table("coverage_proportions.txt",header=TRUE)
cc = read.table("coverage_counts.txt",header=TRUE)
bm = read.table("bam.metadata.fixed",header=FALSE,sep='\t',comment.char='',quote='')
colnames(bm) = c("sname","sid","ge","tech","bam")
big6_snames = read.table("alldata.samples")$V1 # 6 samples that most analyses will be based on 
wgs_snames = read.table("wgs.samples")$V1 # additional WGS samples
# load means from interval-level data
is_be_means = read.table("is_broadexome.bed_all_means.txt",header=TRUE)
is_gc_means = read.table("is_gencode_cds.bed_all_means.txt",header=TRUE)


#### cumulative proportions plots
cp$sname = bm$sname[match(cp$sid,bm$sid)]
cp$ge = bm$ge[match(cp$sid,bm$sid)]
cp$tech = bm$tech[match(cp$sid,bm$sid)]

gte0index = 4 # index of the column containing the "greater than or equal to 0" values
cp$meancov = rep(0.0,dim(cp)[1])
for (n in 1:500) {
  cp$meancov = cp$meancov + n*(cp[,n+gte0index-1]-cp[,n+gte0index])
}

b6 = cp$sname %in% big6_snames
hq = cp$qual == '20_1'
lq = cp$qual == '0_0'
gc = cp$targets == 'gencode_cds'
be = cp$targets == 'broadexome'
cpmat = as.matrix(cp[,4:84])

# mean coverage boxplots
for (target in unique(cp$targets)) {
  for (qual in unique(cp$qual)) {
    png(paste('mean.cov.',target,'.',qual,'.boxplot.png',sep=''),width=1200,height=800)
    subs = cp$qual==qual & cp$targets==target
    boxplot(cp$meancov[subs] ~ paste(cp$ge[subs], cp$tech[subs]),ylim=c(0,120),
            col=c(ecolor,ecolor,gcolor,gcolor),
            main=paste("Mean coverage over",fullname[target],"\n",fullname[qual]))
    categories = paste(cp$ge[subs], cp$tech[subs])
    cat_counts = table(categories)
    mtext(side=1,padj=4,at=1:4,text=paste("n =",cat_counts))
    dev.off()
  } 
}


# cumulative coverage plot for the 6 individuals that have three sequencing datasets
for (target in unique(cp$targets)) {
  for (qual in unique(cp$qual)) {
    png(paste('cumul.cov.',target,'.',qual,'b6.indivs.png',sep=''),width=1200,height=800)
    plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main=paste('Cumulative depth of coverage\nover',fullname[target],"\n",fullname[qual]),
     sub='N = 6 individuals sequenced using all three technologies')
    axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
    for (row in which(b6 & cp$qual==qual & cp$targets==target)) {
      points(cpmat[row,],type='l',lwd=2,
             lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
    }
    legend("bottomleft",legtable$text[1:3],col=legtable$k[1:3],lty=legtable$lty[1:3],lwd=2)
    dev.off()
  }
}

# cumulative coverage plot of HiSeq 2000 vs. X Ten

counts = table(cp$tech[hq & be & cp$ge=='WGS'])
paste(counts, names(counts))
# WGS 2000 vs. X Ten
plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Broad Exome (GATK bundle)',
     sub=paste('N = ',paste(counts, names(counts),collapse=',')))
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(hq & be & cp$ge=='WGS')) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[3:4],col=legtable$k[3:4],lty=legtable$lty[3:4],lwd=2)

# separate them for greater visibility
plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Broad Exome (GATK bundle)',
     sub=paste('N = ',paste(counts, names(counts),collapse=',')))
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(hq & be & cp$ge=='WGS' & cp$tech=='HiSeq 2000')) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[3:4],col=legtable$k[3:4],lty=legtable$lty[3:4],lwd=2)


plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Broad Exome (GATK bundle)',
     sub=paste('N = ',paste(counts, names(counts),collapse=',')))
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(hq & be & cp$ge=='WGS' & cp$tech=='HiSeq X')) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[3:4],col=legtable$k[3:4],lty=legtable$lty[3:4],lwd=2)


# zoom in
plot(NA,NA,xlim=c(0,20),ylim=c(.8,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Broad Exome (GATK bundle)',
     sub='N = ')
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(hq & be & cp$ge=='WGS')) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[3:4],col=legtable$k[3:4],lty=legtable$lty[3:4],lwd=2)


# can add in more analyses from callability project here

#### interval summary stuff
is_be_names = is_be_means[,"interval"]
is_be_chr = sapply(strsplit(is_be_names,":"),"[[",1)
chrbreaks = which(!duplicated(is_be_chr))
chrbreaks = c(chrbreaks,length(is_be_chr))


# mean interval depth by quality on each tech
png('interval.means.ice.by.quality.png',width=1200,height=800)
plot(is_be_means$ICE.0_0,is_be_means$ICE.20_1,pch='.',col=ecolor,
     ylim=c(0,200),xlim=c(0,200),
     xlab='Total depth, any quality',
     ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
     main='Mean Broad Exome interval depth\nin ICE exomes by quality',
     sub=paste('N =',sum(bm$tech=="ICE"),"exomes"))
dev.off()

png('interval.means.agilent.by.quality.png',width=1200,height=800)
plot(is_be_means$Agilent.0_0,is_be_means$Agilent.20_1,pch='.',col=ecolor,
     ylim=c(0,200),xlim=c(0,200),
     xlab='Total depth, any quality',
     ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
     main='Mean Broad Exome interval depth\nin Agilent exomes by quality',
     sub=paste('N =',sum(bm$tech=="Agilent"),"exomes"))
dev.off()

png('interval.means.h2000.by.quality.png',width=1200,height=800)
plot(is_be_means$h2000.0_0,is_be_means$h2000.20_1,pch='.',col=gcolor,
     ylim=c(0,200),xlim=c(0,200),
     xlab='Total depth, any quality',
     ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
     main='Mean Broad Exome interval depth\nin HiSeq 2000 whole genomes by quality',
     sub=paste('N =',sum(bm$tech=="HiSeq 2000"),"genomes"))
dev.off()

png('interval.means.hx.by.quality.png',width=1200,height=800)
plot(is_be_means$hX.0_0,is_be_means$hX.20_1,pch='.',col=gcolor,
     ylim=c(0,200),xlim=c(0,200),
     xlab='Total depth, any quality',
     ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
     main='Mean Broad Exome interval depth\nin HiSeq X Ten whole genomes by quality',
     sub=paste('N =',sum(bm$tech=="HiSeq X"),"genomes"))
dev.off()

# mean interval depth between techs
png('interval.means.ice.vs.agilent.hq.png',width=1200,height=800)
plot(is_be_means$Agilent.20_1,is_be_means$ICE.20_1,pch='.',col=ecolor,
     ylim=c(0,200),xlim=c(0,200),xaxs='i',yaxs='i',
     xlab='Agilent depth',
     ylab='ICE depth',
     main='Mean Broad Exome interval depth\nat BQ ≥ 20 MAPQ ≥ 1\nAgilent vs. ICE',
     sub=paste('N =',sum(bm$tech=="ICE")," ICE exomes vs.",sum(bm$tech=="Agilent")," Agilent exomes"))
m = lm(is_be_means$ICE.20_1 ~ is_be_means$Agilent.20_1)
rho = sqrt(summary(m)$adj.r.squared)
mtext(side=1,text=paste("Pearson's rho =",formatC(rho,digits=2)),col='red',cex=.8)
abline(m,col='red')
dev.off()

png('interval.means.h2.vs.hx.hq.png',width=1200,height=800)
plot(is_be_means$h2000.20_1,is_be_means$hX.20_1,pch='.',col=gcolor,
     ylim=c(0,200),xlim=c(0,200),xaxs='i',yaxs='i',
     xlab='HiSeq 2000 depth',
     ylab='HiSeq X Ten depth',
     main='Mean Broad Exome interval depth\nat BQ ≥ 20 MAPQ ≥ 1\nHiSeq 2000 vs. X Ten',
     sub=paste('N =',sum(bm$tech=="HiSeq 2000")," HiSeq 2000 whole genomes vs.",sum(bm$tech=="HiSeq X")," HiSeq X Ten whole genomes"))
m = lm(is_be_means$hX.20_1 ~ is_be_means$h2000.20_1)
rho = sqrt(summary(m)$adj.r.squared)
mtext(side=1,text=paste("Pearson's rho =",formatC(rho,digits=2)),col='red',cex=.8)
abline(m,col='red')
dev.off()

# which genomic regions account for differences in coverage?
hx_problem_interval_idx = which(is_be_means$hX.20_1 < 20 & is_be_means$hX.20_1 < is_be_means$h2000.20_1)
hx_depth_loss = rep(0.0,dim(is_be_means)[1])
hx_depth_loss[hx_problem_interval_idx] =  -(is_be_means$h2000.20_1[hx_problem_interval_idx] - is_be_means$hX.20_1[hx_problem_interval_idx])/is_be_means$h2000.20_1[hx_problem_interval_idx]
proportion_intervals_problematic = length(hx_problem_interval_idx) / dim(is_be_means)[1]
sub = paste('These "problem" intervals represent ',formatC(proportion_intervals_problematic*100,digits=2),"% of intervals",sep="")
png('depth.loss.regions.xten.2000.png',width=1200,height=800)
par(mar=c(6,6,4,4))
plot(1:dim(is_be_means)[1],hx_depth_loss,pch=20,col=gcolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',
     ylab='% loss of depth in HiSeq X Ten\ncompared to HiSeq 2000',
     xlab='Broad Exome interval',
     main='Intervals where HiSeq X Ten loses high-quality depth compared to HiSeq 2000\nand X Ten has mean depth < 20',
     sub=sub)
abline(v=chrbreaks,col='red')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
axis(side=2,at=-(1:10)/10,labels=paste(-(1:10)*10,"%",sep=""),cex.axis=.8)
dev.off()

hx_problem_interval_idx = which(is_be_means$ICE.20_1 < 20 & is_be_means$ICE.20_1 < is_be_means$Agilent.20_1)
hx_depth_loss = rep(0.0,dim(is_be_means)[1])
hx_depth_loss[hx_problem_interval_idx] =  -(is_be_means$Agilent.20_1[hx_problem_interval_idx] - is_be_means$ICE.20_1[hx_problem_interval_idx])/is_be_means$Agilent.20_1[hx_problem_interval_idx]
proportion_intervals_problematic = length(hx_problem_interval_idx) / dim(is_be_means)[1]
sub = paste('These "problem" intervals represent ',formatC(proportion_intervals_problematic*100,digits=2),"% of intervals",sep="")
png('depth.loss.regions.ice.agilent.png',width=1200,height=800)
par(mar=c(6,6,4,4))
plot(1:dim(is_be_means)[1],hx_depth_loss,pch=20,col=ecolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',
     ylab='% loss of depth in ICE\ncompared to Agilent',
     xlab='Broad Exome interval',
     main='Intervals where ICE loses high-quality depth compared to Agilent\nand ICE has mean depth < 20',
     sub=sub)
abline(v=chrbreaks,col='red')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
axis(side=2,at=-(1:10)/10,labels=paste(-(1:10)*10,"%",sep=""),cex.axis=.8)
dev.off()


#### genotypic concordance

# tables for overall WGS concordance with 2.5M array
gcp_xten_array = read.table("wgs.xten.vs.array.gq30dp10.molt.summary.concordance.proportions",skip=3,header=TRUE)
gcp_2000_array = read.table("wgs.2000.vs.array.gq30dp10.molt.summary.concordance.proportions",skip=3,header=TRUE)

gcp_xten_array_table = acast(data=subset(gcp_xten_array, Comp_Genotype != "Mismatching_Alleles"),
      formula=Eval_Genotype ~ Comp_Genotype,
      value.var="Proportion")
write.table(gcp_xten_array_table,"gcp_xten_array_table.txt",sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)

gcp_2000_array_table = acast(data=subset(gcp_2000_array, Comp_Genotype != "Mismatching_Alleles"),
                             formula=Eval_Genotype ~ Comp_Genotype,
                             value.var="Proportion")
write.table(gcp_2000_array_table,"gcp_2000_array_table.txt",sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
