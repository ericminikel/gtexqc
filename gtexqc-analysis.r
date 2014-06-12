
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

#### read in data
cp = read.table("coverage_proportions.txt",header=TRUE)
cc = read.table("coverage_counts.txt",header=TRUE)
bm = read.table("bam.metadata.fixed",header=FALSE,sep='\t',comment.char='',quote='')
colnames(bm) = c("sname","sid","ge","tech","bam")
big6_snames = read.table("alldata.samples")$V1 # 6 samples that most analyses will be based on 
wgs_snames = read.table("wgs.samples")$V1 # additional WGS samples
is_be = read.table("interval_summary_broadexome.bed.txt",header=TRUE)
#is_gc = read.table("interval_summary_gencode_cds.bed.txt",header=TRUE)
# load means from interval-level data
is_be_means = read.table("is_be_all_means.txt",header=TRUE)


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

boxplot(cp$meancov[lq & gc] ~ paste(cp$ge[lq & gc], cp$tech[lq & gc]),ylim=c(0,120),
        col=c(ecolor,ecolor,gcolor,gcolor),
        main="Mean coverage over Gencode coding intervals",
        sub="No quality filters")
categories = paste(cp$ge[hq & gc], cp$tech[hq & gc])
cat_counts = table(categories)
mtext(side=1,padj=4,at=1:4,text=paste("n =",cat_counts))

boxplot(cp$meancov[hq & gc] ~ paste(cp$ge[hq & gc], cp$tech[hq & gc]),ylim=c(0,120),
        col=c(ecolor,ecolor,gcolor,gcolor),
        main="Mean coverage over Gencode coding intervals",
        sub='BQ ≥ 20 & MAPQ ≥ 1 only')
categories = paste(cp$ge[hq & gc], cp$tech[hq & gc])
cat_counts = table(categories)
mtext(side=1,padj=4,at=1:4,text=paste("n =",cat_counts))



plot(NA,NA,xlim=c(0,80),ylim=c(0,1),xaxs='i',yaxs='i',yaxt='n',
     xlab='Depth',ylab='Proportion of target covered at depth',
     main='Cumulative depth of coverage\nover Broad Exome (GATK bundle)',
     sub='N = 6 individuals sequenced using all three technologies')
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex=.8)
for (row in which(b6 & hq & be)) {
  points(cpmat[row,],type='l',lwd=2,
         lty=getlty(cp$tech[row]),col=getk(cp$ge[row]))
}
legend("bottomleft",legtable$text[1:3],col=legtable$k[1:3],lty=legtable$lty[1:3],lwd=2)

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


## begin deprecated code
is_be_mat = as.matrix(is_be[3:dim(is_be)[2]])

# now that the matrix is separate, add more annotations to the data frame
is_be$sname = bm$sname[match(is_be$sid,bm$sid)]
is_be$ge    = bm$ge[match(is_be$sid,bm$sid)]
is_be$tech  = bm$tech[match(is_be$sid,bm$sid)]
is_be_meta = is_be[,c("sid","qual","sname","ge","tech")]

# get a list of sample ids for which both high and low quality counts are available
has_both = unique(is_be$sid[duplicated(is_be$sid)]) 

# filters
hb = is_be_meta$sid %in% has_both
b6 = is_be_meta$sname %in% big6_snames
hq = is_be_meta$qual == '20_1'
lq = is_be_meta$qual == '0_0'
wg = is_be_meta$ge == "WGS"
we = is_be_meta$ge == "WES"
ag = is_be_meta$tech == "Agilent"
ic = is_be_meta$tech == "ICE"
h2 = is_be_meta$tech == "HiSeq 2000"
hx = is_be_meta$tech == "HiSeq X"

# first plot average depth histograms
# separately
for (tech in unique(is_be_meta$tech)) {
  techstring = gsub(" ","",tolower(tech))
  png(paste('interval.avdepth.',techstring,".be.hq.png",sep=""),width=1200,height=800)
  is_be_subs  = as.matrix(is_be_mat[b6 & hq & hb & is_be_meta$tech==tech,])
  is_be_subs_means = colMeans(is_be_subs)
  plot(1:dim(is_be_mat)[2],is_be_subs_means,pch='.',col=ecolor,
       xaxs='i',,yaxs='i',xaxt='n',xlab='',ylim=c(0,200),
       ylab='Mean depth',
       main=paste('Mean depth by Broad Exome interval\n',tech))
  axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
  abline(v=chrbreaks,col='red')
  abline(h=0,col='black',lwd=2)
  dev.off()
}
# is_be_ag  = as.matrix(is_be_mat[b6 & hq & hb & we & ag,])
# is_be_ag_means = colMeans(is_be_ag)
# plot(1:dim(is_be_mat)[2],is_be_ag_means,pch='.',col=ecolor,
#      xaxs='i',,yaxs='i',xaxt='n',xlab='',ylim=c(0,200),
#      ylab='Mean depth',
#      main='Mean depth by Broad Exome interval\nAgilent Exomes')
# axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
# abline(v=chrbreaks,col='red')
# abline(h=0,col='black',lwd=2)
# 
# is_be_ic  = as.matrix(is_be_mat[b6 & hq & hb & we & ic,])
# is_be_ic_means = colMeans(is_be_ic)
# plot(1:dim(is_be_mat)[2],is_be_ic_means,pch='.',col=ecolor,
#      xaxs='i',,yaxs='i',xaxt='n',xlab='',ylim=c(0,200),
#      ylab='Mean depth',
#      main='Mean depth by Broad Exome interval\nICE Exomes')
# axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
# abline(v=chrbreaks,col='red')
# abline(h=0,col='black',lwd=2)

# also try combined plot pointing in opposite directions
plot(1:dim(is_be_mat)[2],is_be_ag_means,pch='.',col=ecolor,
     xaxs='i',,yaxs='i',xaxt='n',xlab='',ylim=c(-200,200),yaxt='n',
     ylab='Mean depth',
     main='Mean depth by Broad Exome interval\nAgilent Exomes')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
abline(v=chrbreaks,col='red')
abline(h=0,col='black',lwd=2)
points(1:dim(is_be_mat)[2],-is_be_ic_means,pch='.',col=ecolor)
axis(side=2,at=c(-200,-100,0,100,200),labels=c(200,100,0,100,200))
mtext(side=2,at=c(-150,150),padj=-2,text=c("ICE","Agilent"))


# find the ratio of high-quality to all depth by interval,
# in only the samples with both
# first do this for WES Agilent
is_be_20_1  = as.matrix(is_be_mat[hq & hb & we & ag,])
is_be_0_0   = as.matrix(is_be_mat[lq & hb & we & ag,])
is_be_ratio = is_be_20_1 / is_be_0_0
is_be_ratio_interval_mean = colMeans(is_be_ratio)
plot(1:dim(is_be_mat)[2],is_be_ratio_interval_mean,pch='.',col=ecolor,
     xaxs='i',,yaxs='i',xaxt='n',xlab='',
     ylab='ratio',
     main='Ratio of MAPQ ≥ 20 & BQ ≥ 1 coverage\nto all coverage, by interval\nAgilent Exomes')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
abline(v=chrbreaks,col='red')

is_be_20_1  = as.matrix(is_be_mat[hq & hb & we & ic,])
is_be_0_0   = as.matrix(is_be_mat[lq & hb & we & ic,])
is_be_ratio = is_be_20_1 / is_be_0_0
is_be_ratio_interval_mean = colMeans(is_be_ratio)
plot(1:dim(is_be_mat)[2],is_be_ratio_interval_mean,pch='.',col=ecolor,
     xaxs='i',,yaxs='i',xaxt='n',xlab='',
     ylab='ratio',
     main='Ratio of MAPQ ≥ 20 & BQ ≥ 1 coverage\nto all coverage, by interval\nICE Exomes')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
abline(v=chrbreaks,col='red')



# then for WGS
is_be_20_1  = as.matrix(is_be_mat[hq & hb & wg,])
is_be_0_0   = as.matrix(is_be_mat[lq & hb & wg,])
is_be_ratio = is_be_20_1 / is_be_0_0
is_be_ratio_interval_mean = colMeans(is_be_ratio)


plot(1:dim(is_be_mat)[2],is_be_ratio_interval_mean,pch='.',col=ecolor,
     xaxs='i',,yaxs='i',xaxt='n',xlab='',
     ylab='ratio',
     main='Ratio of MAPQ ≥ 20 & BQ ≥ 1 coverage\nto all coverage, by interval')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
abline(v=chrbreaks,col='red')



has_both = unique(is_be$sid[duplicated(is_be$sid)]) 
hb = is_be_meta$sid %in% has_both
b6 = is_be_meta$sname %in% big6_snames
hq = is_be_meta$qual == '20_1'
is_be_ice = is_be_mat[hb & b6 & hq & is_be_meta$tech=="ICE",]
is_be_ag  = is_be_mat[hb & b6 & hq & is_be_meta$tech=="Agilent",]
is_be_ice_means = colMeans(is_be_ice)
is_be_ag_means = colMeans(is_be_ag)
# plot where each point is the mean of all samples for 1 interval
plot(is_be_ice_means,is_be_ag_means,pch='.',col=ecolor,
     xlim=c(0,80),ylim=c(0,80))
abline(a=0,b=1,col='red')
cor.test(is_be_ice_means,is_be_ag_means) # r = .19, p < 2e-16
# plot where each point is 1 sample 1 interval?


# end deprecated code

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

