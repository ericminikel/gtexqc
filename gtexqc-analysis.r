
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


#### cumulative proportions plots
cp$sname = bm$sname[match(cp$sid,bm$sid)]
cp$ge = bm$ge[match(cp$sid,bm$sid)]
cp$tech = bm$tech[match(cp$sid,bm$sid)]

gte0index = 4 # index of the column containing the "greater than or equal to 0" values
cp$meancov = rep(0.0,dim(cp)[1])
for (n in 1:500) {
  cp$meancov = cp$meancov + n*(cp[,n+gte0index-1]-cp[,n+gte0index])
}

boxplot(cp$meancov ~ paste(cp$ge, cp$tech),ylim=c(0,120))

b6 = cp$sname %in% big6_snames
hq = cp$qual == '20_1'
gc = cp$targets == 'gencode_cds'
be = cp$targets == 'broadexome'
cpmat = as.matrix(cp[4:84])
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

# can add in more analyses from callability project here

#### interval summary stuff

is_be_names = colnames(is_be[3:dim(is_be)[2]]) # get interval names
is_be_chr = sapply(strsplit(substr(is_be_names,2,nchar(is_be_names)),"\\."),"[[",1) # extract chromosome
chrbreaks = which(!duplicated(is_be_chr)) # find the first unique instance of each chromosome
chrbreaks = c(chrbreaks,length(is_be_chr)) # add a final break for end of last chromosome

is_be_mat = as.matrix(is_be[3:dim(is_be)[2]])

# now that the matrix is separate, add more annotations to the data frame
is_be$sname = bm$sname[match(is_be$sid,bm$sid)]
is_be$ge    = bm$ge[match(is_be$sid,bm$sid)]
is_be$tech  = bm$tech[match(is_be$sid,bm$sid)]

is_be_20_1 = colMeans(is_be_mat[is_be$qual=="20_1",])
is_be_0_0 =  colMeans(is_be_mat[is_be$qual=="0_0",])

# get a list of sample ids for which both high and low quality counts are available
has_both = unique(is_be$sid[duplicated(is_be$sid)]) 

# find the ratio of high-quality to all depth by interval,
# in only the samples with both
is_be_20_1  = as.matrix(is_be[is_be$qual=="20_1" & is_be$sid %in% has_both,3:dim(is_be)[2]])
is_be_0_0   = as.matrix(is_be[is_be$qual=="0_0"  & is_be$sid %in% has_both,3:dim(is_be)[2]])
is_be_ratio = is_be_20_1 / is_be_0_0
is_be_ratio_interval_mean = colMeans(is_be_ratio)

plot(3:dim(is_be)[2],is_be_ratio_interval_mean,pch='.',col='blue',
     xaxs='i',,yaxs='i',xaxt='n',yaxt='n',xlab='',
     ylab='ratio',
     main='Ratio of MAPQ ≥ 20 & BQ ≥ 1 coverage\nto all coverage, by interval')
axis(side=1,at=midpoints(chrbreaks),labels=unique(is_be_chr),lty=0,cex.axis=.8)
abline(v=chrbreaks,col='red')
