
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

getk = function(ge_or_tech) {
  if (ge_or_tech %in% c("WES","ICE","Agilent")) {
    return (ecolor)
  } else if (ge_or_tech %in% c("WGS","HiSeq 2000","HiSeq X")) {
    return (gcolor)
  } else {
    stop("invalid specification - see getk() function")
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
is_means_meta=data.frame(tech=c("","ICE","ICE","Agilent","Agilent","HiSeq 2000","HiSeq 2000","HiSeq X","HiSeq X"),
                            qual=c("","20_1","0_0","20_1","0_0","20_1","0_0","20_1","0_0"))

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
    png(paste('cumul.cov.',target,'.',qual,'.b6.indivs.png',sep=''),width=1200,height=800)
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

median(cp$gte_10[cp$qual=="20_1" & cp$tech=="HiSeq X" & cp$targets=="gencode_cds" & b6])
median(cp$gte_10[cp$qual=="20_1" & cp$tech=="ICE" & cp$targets=="gencode_cds" & b6])
median(cp$gte_10[cp$qual=="20_1" & cp$tech=="Agilent" & cp$targets=="gencode_cds" & b6])

median(cp$gte_20[cp$qual=="20_1" & cp$tech=="HiSeq X" & cp$targets=="gencode_cds" & b6])
median(cp$gte_20[cp$qual=="20_1" & cp$tech=="ICE" & cp$targets=="gencode_cds" & b6])
median(cp$gte_20[cp$qual=="20_1" & cp$tech=="Agilent" & cp$targets=="gencode_cds" & b6])


# cumulative coverage plot of HiSeq 2000 vs. X Ten
for (target in unique(cp$targets)) {
  for (qual in unique(cp$qual)) {
    png(paste('cumul.cov.',target,'.',qual,'.wgs.png',sep=''),width=1200,height=800)
    counts = table(cp$tech[cp$qual==qual & cp$targets==target & cp$ge=='WGS'])
    paste(counts, names(counts))
    # WGS 2000 vs. X Ten
    plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
         xlab='Depth',ylab='Proportion of target covered at depth',
         main=paste('Cumulative depth of coverage\nover ',fullname[target],'\n',fullname[qual],sep=''),
         sub=paste('N = ',paste(counts, names(counts),collapse=',')))
    axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
    for (row in which(cp$qual==qual & cp$targets==target & cp$ge=='WGS')) {
      points(cpmat[row,],type='l',lwd=2,
             lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
    }
    legend("bottomleft",legtable$text[3:4],col=legtable$k[3:4],lty=legtable$lty[3:4],lwd=2)
    dev.off() 
  }
}

# separate them for greater visibility
qual = "20_1"
ge = "WGS"
target = "gencode_cds"

png('cumul.cov.h2000.gencode_cds.20_1.png',width=1200,height=800)
counts = table(cp$tech[cp$tech=='HiSeq 2000' & cp$qual==qual & cp$targets==target & cp$ge==ge])
plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Gencode CDS intervals',
     sub=paste('N = ',paste(counts, names(counts),collapse=',')))
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(cp$tech=='HiSeq 2000' & cp$qual==qual & cp$targets==target & cp$ge==ge)) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[4],col=legtable$k[4],lty=legtable$lty[4],lwd=2)
dev.off()

png('cumul.cov.hx.gencode_cds.20_1.png',width=1200,height=800)
counts = table(cp$tech[cp$tech=='HiSeq X' & cp$qual==qual & cp$targets==target & cp$ge==ge])
plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Gencode CDS intervals',
     sub=paste('N = ',paste(counts, names(counts),collapse=',')))
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(cp$tech=='HiSeq X' & cp$qual==qual & cp$targets==target & cp$ge==ge)) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[3],col=legtable$k[3],lty=legtable$lty[3],lwd=2)
dev.off()

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
be_chrbreaks = which(!duplicated(is_be_chr))
be_chrbreaks = c(be_chrbreaks,length(is_be_chr))

is_gc_names = is_gc_means[,"interval"]
is_gc_chr = sapply(strsplit(is_gc_names,":"),"[[",1)
gc_chrbreaks = which(!duplicated(is_gc_chr))
gc_chrbreaks = c(gc_chrbreaks,length(is_gc_chr))


# mean interval depth by quality on each tech
for (target in unique(cp$targets)) {
  if (target == "broadexome") {
    is_table = is_be_means
  } else if (target == "gencode_cds") {
    is_table = is_gc_means
  }
  png(paste('interval.means.ice.by.quality.',target,'.png',sep=''),width=1200,height=800)
  plot(is_table$ICE.0_0,is_table$ICE.20_1,pch='.',col=ecolor,
       ylim=c(0,500),xlim=c(0,500),
       xlab='Total depth, any quality',
       ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
       main=paste('Mean ',fullname[target],' interval depth\nin ICE exomes by quality',sep=''),
       sub=paste('N =',sum(bm$tech=="ICE"),"exomes"))
  dev.off()
  
  png(paste('interval.means.agilent.by.quality.',target,'.png',sep=''),width=1200,height=800)
  plot(is_table$Agilent.0_0,is_table$Agilent.20_1,pch='.',col=ecolor,
       ylim=c(0,500),xlim=c(0,500),
       xlab='Total depth, any quality',
       ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
       main=paste('Mean ',fullname[target],' interval depth\nin Agilent exomes by quality',sep=''),
       sub=paste('N =',sum(bm$tech=="Agilent"),"exomes"))
  dev.off()
  
  png(paste('interval.means.h2000.by.quality.',target,'.png',sep=''),width=1200,height=800)
  plot(is_table$h2000.0_0,is_table$h2000.20_1,pch='.',col=gcolor,
       ylim=c(0,200),xlim=c(0,200),
       xlab='Total depth, any quality',
       ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
       main=paste('Mean ',fullname[target],' interval depth\nin HiSeq 2000 genomes by quality',sep=''),
       sub=paste('N =',sum(bm$tech=="HiSeq 2000"),"genomes"))
  dev.off()
  
  png(paste('interval.means.hx.by.quality.',target,'.png',sep=''),width=1200,height=800)
  plot(is_table$hX.0_0,is_table$hX.20_1,pch='.',col=gcolor,
       ylim=c(0,200),xlim=c(0,200),
       xlab='Total depth, any quality',
       ylab='Depth with BQ ≥ 20, MAPQ ≥ 1',
       main=paste('Mean ',fullname[target],' interval depth\nin HiSeq X Ten genomes by quality',sep=''),
       sub=paste('N =',sum(bm$tech=="HiSeq X"),"genomes"))
  dev.off()
}


# plot across the genome
for (tech in unique(bm$tech)) {
  qual = "20_1"
  png(paste('chrom.coord.mean.cov.',tech,'.20_1.png',sep=''),width=1200,height=800)
      usecol = which(is_means_meta$tech==tech & is_means_meta$qual==qual)
      par(mar=c(6,6,4,4))
      plot(1:dim(is_gc_means)[1],is_gc_means[,usecol],pch='.',col=getk(tech),
           xaxt='n',xaxs='i',yaxs='i',yaxt='n',ylim=c(0,200),
           ylab='Mean depth',
           xlab='Interval',
           main=paste('Mean depth by Gencode CDS interval\n',fullname[tech]," ",fullname[qual],sep=''))
      abline(v=gc_chrbreaks,col='red')
      axis(side=1,at=midpoints(gc_chrbreaks),labels=unique(is_gc_chr),lty=0,cex.axis=.8)
      axis(side=2,at=c(50,100,150,200),labels=c(50,100,150,200),cex.axis=.8)
      dev.off()
}


# mean interval depth between techs
png('interval.means.ice.vs.agilent.20_1.broadexome.png',width=1200,height=800)
plot(is_be_means$Agilent.20_1,is_be_means$ICE.20_1,pch='.',col=ecolor,
     ylim=c(0,200),xlim=c(0,200),xaxs='i',yaxs='i',
     xlab='Agilent depth',
     ylab='ICE depth',
     main='Mean BroadExome interval depth\nat BQ ≥ 20 MAPQ ≥ 1\nAgilent vs. ICE',
     sub=paste('N =',sum(bm$tech=="ICE")," ICE exomes vs.",sum(bm$tech=="Agilent")," Agilent exomes"))
m = lm(is_be_means$ICE.20_1 ~ is_be_means$Agilent.20_1)
rho = sqrt(summary(m)$adj.r.squared)
mtext(side=1,text=paste("Pearson's rho =",formatC(rho,digits=2)),col='red',cex=.8)
abline(a=0,b=1,col='#000000')
dev.off()

png('interval.means.ice.vs.agilent.20_1.gencode_cds.png',width=1200,height=800)
plot(is_gc_means$Agilent.20_1,is_gc_means$ICE.20_1,pch='.',col=ecolor,
     ylim=c(0,200),xlim=c(0,200),xaxs='i',yaxs='i',
     xlab='Agilent depth',
     ylab='ICE depth',
     main='Mean Gencode CDS interval depth\nat BQ ≥ 20 MAPQ ≥ 1\nAgilent vs. ICE',
     sub=paste('N =',sum(bm$tech=="ICE")," ICE exomes vs.",sum(bm$tech=="Agilent")," Agilent exomes"))
m = lm(is_gc_means$ICE.20_1 ~ is_gc_means$Agilent.20_1)
rho = sqrt(summary(m)$adj.r.squared)
mtext(side=1,text=paste("Pearson's rho =",formatC(rho,digits=2)),col='red',cex=.8)
abline(a=0,b=1,col='#000000')
dev.off()

png('interval.means.h2.vs.hx.20_1.broadexome.png',width=1200,height=800)
plot(is_be_means$h2000.20_1,is_be_means$hX.20_1,pch='.',col=gcolor,
     ylim=c(0,200),xlim=c(0,200),xaxs='i',yaxs='i',
     xlab='HiSeq 2000 depth',
     ylab='HiSeq X Ten depth',
     main='Mean Broad Exome interval depth\nat BQ ≥ 20 MAPQ ≥ 1\nHiSeq 2000 vs. X Ten',
     sub=paste('N =',sum(bm$tech=="HiSeq 2000")," HiSeq 2000 whole genomes vs.",sum(bm$tech=="HiSeq X")," HiSeq X Ten whole genomes"))
m = lm(is_be_means$hX.20_1 ~ is_be_means$h2000.20_1)
rho = sqrt(summary(m)$adj.r.squared)
mtext(side=1,text=paste("Pearson's rho =",formatC(rho,digits=2)),col='red',cex=.8)
abline(a=0,b=1,col='#000000')
dev.off()

png('interval.means.h2.vs.hx.20_1.gencode_cds.png',width=1200,height=800)
plot(is_gc_means$h2000.20_1,is_gc_means$hX.20_1,pch='.',col=gcolor,
     ylim=c(0,200),xlim=c(0,200),xaxs='i',yaxs='i',
     xlab='HiSeq 2000 depth',
     ylab='HiSeq X Ten depth',
     main='Mean Gencode CDS interval depth\nat BQ ≥ 20 MAPQ ≥ 1\nHiSeq 2000 vs. X Ten',
     sub=paste('N =',sum(bm$tech=="HiSeq 2000")," HiSeq 2000 whole genomes vs.",sum(bm$tech=="HiSeq X")," HiSeq X Ten whole genomes"))
m = lm(is_gc_means$hX.20_1 ~ is_gc_means$h2000.20_1)
rho = sqrt(summary(m)$adj.r.squared)
mtext(side=1,text=paste("Pearson's rho =",formatC(rho,digits=2)),col='red',cex=.8)
abline(a=0,b=1,col='#000000')
dev.off()



# which genomic regions account for differences in coverage?
hx_problem_interval_idx = which(is_gc_means$ICE.20_1 < 20 & is_gc_means$ICE.20_1 < is_gc_means$Agilent.20_1)
hx_depth_loss = rep(0.0,dim(is_gc_means)[1])
hx_depth_loss[hx_problem_interval_idx] =  -(is_gc_means$Agilent.20_1[hx_problem_interval_idx] - is_gc_means$ICE.20_1[hx_problem_interval_idx])/is_gc_means$Agilent.20_1[hx_problem_interval_idx]
proportion_intervals_problematic = length(hx_problem_interval_idx) / dim(is_gc_means)[1]
sub = paste('These "problem" intervals represent ',formatC(proportion_intervals_problematic*100,digits=2),"% of intervals",sep="")
png('depth.loss.regions.ice.agilent.20_1.png',width=1200,height=800)
par(mar=c(6,6,4,4))
plot(1:dim(is_gc_means)[1],hx_depth_loss,pch=20,col=ecolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',
     ylab='% loss of depth in ICE\ncompared to Agilent',
     xlab='Gencode CDS interval',
     main='Intervals where ICE loses high-quality depth compared to Agilent\nand ICE has mean depth < 20',
     sub=sub)
abline(v=gc_chrbreaks,col='red')
axis(side=1,at=midpoints(gc_chrbreaks),labels=unique(is_gc_chr),lty=0,cex.axis=.8)
axis(side=2,at=-(1:10)/10,labels=paste(-(1:10)*10,"%",sep=""),cex.axis=.8)
dev.off()


hx_problem_interval_idx = which(is_gc_means$hX.20_1 < 20 & is_gc_means$hX.20_1 < is_gc_means$h2000.20_1)
hx_depth_loss = rep(0.0,dim(is_gc_means)[1])
hx_depth_loss[hx_problem_interval_idx] =  -(is_gc_means$h2000.20_1[hx_problem_interval_idx] - is_gc_means$hX.20_1[hx_problem_interval_idx])/is_gc_means$h2000.20_1[hx_problem_interval_idx]
proportion_intervals_problematic = length(hx_problem_interval_idx) / dim(is_gc_means)[1]
sub = paste('These "problem" intervals represent ',formatC(proportion_intervals_problematic*100,digits=2),"% of intervals",sep="")
png('depth.loss.regions.xten.2000.20_1.png',width=1200,height=800)
par(mar=c(6,6,4,4))
plot(1:dim(is_gc_means)[1],hx_depth_loss,pch=20,col=gcolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',
     ylab='% loss of depth in HiSeq X Ten\ncompared to HiSeq 2000',
     xlab='Gencode CDS interval',
     main='Intervals where HiSeq X Ten loses high-quality depth compared to HiSeq 2000\nand X Ten has mean depth < 20',
     sub=sub)
abline(v=gc_chrbreaks,col='red')
axis(side=1,at=midpoints(gc_chrbreaks),labels=unique(is_gc_chr),lty=0,cex.axis=.8)
axis(side=2,at=-(1:10)/10,labels=paste(-(1:10)*10,"%",sep=""),cex.axis=.8)
dev.off()


# apply adjustment of +11% to X chromosome for the different male-female ratio in the X Ten individuals
is_gc_means_hX_20_1_adjusted = is_gc_means$hX.20_1
is_gc_means_hX_20_1_adjusted[is_gc_chr=='X'] = is_gc_means_hX_20_1_adjusted[is_gc_chr=='X']*(1.42/1.28)
hx_problem_interval_idx = which(is_gc_means_hX_20_1_adjusted< 20 & is_gc_means_hX_20_1_adjusted < is_gc_means$h2000.20_1)
hx_depth_loss = rep(0.0,dim(is_gc_means)[1])
hx_depth_loss[hx_problem_interval_idx] =  -(is_gc_means$h2000.20_1[hx_problem_interval_idx] - is_gc_means_hX_20_1_adjusted[hx_problem_interval_idx])/is_gc_means$h2000.20_1[hx_problem_interval_idx]
proportion_intervals_problematic = length(hx_problem_interval_idx) / dim(is_gc_means)[1]
sub = paste('These "problem" intervals represent ',formatC(proportion_intervals_problematic*100,digits=2),"% of intervals",sep="")
png('depth.loss.regions.xten.2000.20_1.sex_adjusted.png',width=1200,height=800)
par(mar=c(6,6,4,4))
plot(1:dim(is_gc_means)[1],hx_depth_loss,pch=20,col=gcolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',
     ylab='% loss of depth in HiSeq X Ten\ncompared to HiSeq 2000',
     xlab='Gencode CDS interval',
     main='Intervals where HiSeq X Ten loses high-quality depth compared to HiSeq 2000\nand X Ten has mean depth < 20\nAdjusted for sex ratio',
     sub=sub)
abline(v=gc_chrbreaks,col='red')
axis(side=1,at=midpoints(gc_chrbreaks),labels=unique(is_gc_chr),lty=0,cex.axis=.8)
axis(side=2,at=-(1:10)/10,labels=paste(-(1:10)*10,"%",sep=""),cex.axis=.8)
dev.off()

#### genotypic concordance

# summary
genomewide_gc = read.table("all.concordance.summary",header=TRUE,row.names=1)
genomewide_gc = genomewide_gc[!(grepl("exchip",rownames(genomewide_gc))),-1]
interval_gc = read.table("interval.concordance.summary",header=TRUE,row.names=1)
broadexome_gc = interval_gc[grepl("broadexome",rownames(interval_gc)),-1]
gencode_gc  = interval_gc[grepl("gencode",rownames(interval_gc)),-1]
colnames(genomewide_gc) = paste('genomewide.',colnames(genomewide_gc),sep='')
colnames(broadexome_gc) = paste('broadexome.',colnames(broadexome_gc),sep='')
colnames(gencode_gc) = paste('gencode.',colnames(gencode_gc),sep='')
master_gc_table = rbind(t(genomewide_gc),t(broadexome_gc),t(gencode_gc))

write.table(master_gc_table,"master_gc_table.txt",row.names=TRUE,col.names=TRUE,quote=FALSE)

# detailed
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


# is higher no call rate at hom ref sites in X Ten due to GQ and DP criteria???
gcp_xten_array_noqc = read.table("gc.wgs.xten.vs.array.gencode_cds.bed.molt2.concordance.proportions",skip=3,header=TRUE)
gcp_xten_array_noqc_table = acast(data=subset(gcp_xten_array_noqc, Comp_Genotype != "Mismatching_Alleles"),
                             formula=Eval_Genotype ~ Comp_Genotype,
                             value.var="Proportion")
write.table(gcp_xten_array_noqc_table,"gcp_xten_array_noqc_table.txt",sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)

