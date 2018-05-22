## 4) Processing Scan Results

Processing the scan results is done in R.  Assuming you have installed qtl2 on your local computer, this part can be done on your local computer.  Much of the output from this part will be plots/graphs and data structures which you will want to examine.  This is not easily done on the linux based servers unless you are using a virtual desktop like VNC, which can be a bit cumbersome.  It is also likely that working interactively will be better, as results from one step will often inform/alter the subsequent steps. Therefore, the recommended approach is to do this locally and interactively in RStudio.  All further instruction will assume you are working locally/interactively in RStudio.

##### 1) Transfer the output from the scans to a directory your local computer.  The simplest approach is using a file transfer program like Filezilla, although other methods like Rsync can be utilized.  pr.Rdata is usually quite a big file and may take some time to transfer.

What do these R objects represent? Unfortunately, these objects are not well documented and you will have play around with them to familiarize yourself with what's in them and how to access the parts that you want.  

There is a good [tutorial on QTL mapping](https://smcclatchy.github.io/mapping/aio/) by Susan McClatchy, which also examines the structure of some of these data structures. Working through this tutorial would give you a nice understanding of QTL mapping and the data structures in qtl2.  

There is also some info on the data structures on [Karl's developer guide page ](http://kbroman.org/qtl2/assets/vignettes/developer_guide.html#qtl_data_structure), but the tutorial mentioned above gives a more comprehensive picture.    

The code below from the tutorial examines some of the data structures and then demonstrates the actual processing of the data to view the LOD plots, find peaks, 95% confidence interval(CI) for peaks, genes within the 95% CI, allele probabilites at the peaks, etc.
```
# set the working directory
setwd("~/path/to/Rqtl2_tutorial")
# load the qtl2 library
library(qtl2)
# reference a file with some convenience functions
# put this in your working directory
# you can load the reference file and look at the functions
# *** NOTE that these convenience functions are not part of the Rqtl2, or any *** #
# *** other published package.  They were written by a member of the Sethupathy *** #
# *** lab to perform repetitive tasks, i.e. to avoid having redundant code ***
# *** in your R file.  Feel free to copy the code and modify it to fit your needs. *** #
source("processQtl2Functions.R")

# load the output files from the scan
load("out.Rdata")
load("do.Rdata")
load("apr.Rdata")
load("pr.Rdata")
load("perms.Rdata")

################################################## 
      #Quick Intro to Data Structures in R/qtl2 
#################################################

###### out.Rdata #############
# out.Rdata contains the LOD scores for each marker/phenotype
colnames(out) # gives us phenotype names
# [1] "fc"    "ob"    "ppi6"  "ppi6F"

# this will give us the first 10 rows, showing the
# marker names in far left col, then LOD scores
# for each pheno at that marker
head(out, n=10)
#                 fc       ob     ppi6    ppi6F
# JAX00000010  1.252868 1.164398 1.024076 1.423277
# JAX00000010r 1.252868 1.164398 1.024076 1.423277
# UNCHS000003  1.252868 1.164399 1.024076 1.423277
# JAX00240606  1.252868 1.164399 1.024077 1.423278
# JAX00240613  1.252867 1.164400 1.024077 1.423279
# JAX00240636  1.252850 1.164373 1.024086 1.423301
# UNCHS000005  1.251154 1.148089 1.026729 1.435749
# UNCHS000006  1.251036 1.146919 1.026973 1.436690
# JAX00240649  1.250465 1.141304 1.028319 1.441402
# UNCHS000007  1.247499 1.109692 1.039076 1.470550

########### do.Rdata ###############
# the do object contains information about our cross
# names of attributes in do
names(do)
#[1] "crosstype"    "geno"         "gmap"         "pmap"         "pheno"        "covar"       
#[7] "founder_geno" "is_x_chr"     "is_female"    "cross_info"   "alleles" 

# do$geno is a 3 dimensional array much like an excel workbook
# each "page" contains the data for one chromosome
# the page is a 2D a matrix with the samples(mice) as the rows and
# markers as the column
# list genotypes for first 5 mice at the first 3 markers on chromosome `1`
do$geno$`1`[1:5,1:3]
#         JAX00000010 JAX00000010r UNCHS000003
# 1            3            3           2
# 13           2            2           3
# 25           3            3           1
# 37           1            1           1
# 2            1            1           3

# do$founder contains the founder gentoypes
do$founder_geno$`1`[1:5,1:3]
#         JAX00000010 JAX00000010r UNCHS000003
# A           1            1           1
# B           1            1           1
# C           3            3           1
# D           3            3           1
# E           1            1           1

#do gmap and pmap contain the genetic and physical location of the
# markers on each chromosome in cM or Mbp, respectively
names(do$gmap)
# [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X"
# list t he locations for the first 5 markers on chr`1`
do$gmap$`1`[1:5]
#JAX00000010 JAX00000010r  UNCHS000003  JAX00240606  JAX00240613 
#0.01227467   0.01227467   0.01773417   0.02381576   0.02931388 

#look at the location in Mbp for the first 3 markers on chr 1
head(do$pmap$`1`,n=3)
#JAX00000010 JAX00000010r  UNCHS000003 
#3.135418     3.135418     3.195649 

# do$pheno is a 2D array of values for phenotypes of each mouse
# rows are mice and columns are phenotype
do$pheno[1:5,] # values for first five mice
#       fc       ob     ppi6    ppi6F
# 1 3.049843 1.000000 1.764206 1.764206
# 2 1.948987 1.300694 1.671304 1.671304
# 3 1.453318 1.389718 1.957862 1.957862
# 4 2.807866 1.405131 1.847511 1.847511
# 5 3.305116 1.360058 1.806612 1.806612
do$covar[1:5,] # covar values for first 5 mice
# sex ngen
# 1   M   25
# 2   M   25
# 3   M   25
# 4   M   25
# 5   M   25

########### pr.Rdata and apr.Rdata #############
# contain the genotype and allele probabilities
# Each 3d array of probabilities is arranged as individuals × genotypes × positions
# look at the names in the 3 dimensions
dimnames(pr$`1`)
# look at the probabilities for the first five mice at marker "UNCHS000005"
pr$`1`[1:5,,"UNCHS000005"]
# AA           AB           BB           AC           BC           CC           AD           BD
# 1  4.834295e-19 9.649268e-19 4.814974e-19 4.256663e-14 4.156218e-14 1.400966e-14 2.574670e-11 2.574569e-11
# 13 1.008685e-19 2.005398e-19 9.967132e-20 1.538924e-17 1.416058e-17 1.647028e-18 1.538924e-17 1.416058e-17
# 25 4.271903e-16 8.111444e-16 3.853032e-16 9.365682e-09 9.329222e-09 1.027755e-03 9.365652e-09 9.329199e-09
# 37 5.084162e-07 1.423987e-03 9.971499e-01 1.991953e-11 2.789730e-08 8.609048e-16 1.991953e-11 2.789730e-08
# 2  8.915764e-16 1.533374e-15 6.623907e-16 5.010867e-16 4.952915e-16 1.462735e-14 5.010867e-16 4.952915e-16
# CD           DD           AE           BE           CE           DE           EE           AF
# 1  1.361106e-11 1.359705e-11 9.668589e-19 9.649268e-19 4.256663e-14 2.574670e-11 4.834295e-19 9.374150e-09
# 13 3.294056e-18 1.647028e-18 2.017370e-19 2.005398e-19 1.538924e-17 1.538924e-17 1.008685e-19 1.955512e-09
# 25 9.979451e-01 1.026565e-03 8.543805e-16 8.111444e-16 9.365682e-09 9.365652e-09 4.271903e-16 3.318690e-19
# 37 1.721809e-15 8.609046e-16 1.016832e-06 1.423987e-03 1.991953e-11 1.991953e-11 5.084162e-07 1.100413e-13
# 2  2.925470e-14 1.462735e-14 1.783153e-15 1.533374e-15 5.010867e-16 5.010867e-16 8.915764e-16 1.005266e-20
# BF           CF           DF           EF           FF           AG           BG           CG
# 1  9.337653e-09 1.029144e-03 9.989705e-01 9.374150e-09 2.721789e-11 3.868384e-15 3.860503e-15 1.423260e-10
# 13 1.932896e-09 1.209211e-07 1.209211e-07 1.955512e-09 3.630387e-12 5.063970e-16 5.033243e-16 2.480034e-14
# 25 3.260676e-19 1.361159e-11 1.361156e-11 3.318690e-19 3.236174e-19 2.043264e-18 1.955856e-18 3.114218e-11
# 37 1.541132e-10 2.158113e-18 2.158113e-18 1.100413e-13 6.076848e-18 6.278161e-13 8.792588e-10 1.232598e-17
# 2  9.544189e-21 4.395385e-19 4.395385e-19 1.005266e-20 3.375324e-21 2.456459e-20 2.166592e-20 4.516445e-19
# DG           EG           FG           GG           AH           BH           CH           DH
# 1  1.063897e-07 3.868384e-15 3.656747e-11 6.270161e-18 2.267936e-17 2.204039e-17 3.258740e-15 1.819925e-12
# 13 2.480034e-14 5.063970e-16 1.868272e-12 3.971825e-19 2.577167e-11 2.576866e-11 6.818433e-12 6.818433e-12
# 25 3.114210e-11 2.043264e-18 1.467754e-18 2.755474e-16 4.532433e-15 4.514784e-15 2.419061e-07 2.419055e-07
# 37 1.232598e-17 6.278161e-13 6.785523e-17 1.937637e-16 8.225917e-12 1.152037e-08 7.086247e-16 7.086247e-16
# 2  4.516445e-19 2.456459e-20 7.804579e-21 6.921859e-21 3.909235e-09 3.864025e-09 2.416251e-07 2.416251e-07
# EH           FH           GH           HH
# 1  2.267936e-17 2.421785e-07 7.138709e-14 8.821286e-19
# 13 2.577167e-11 9.999996e-01 1.064136e-07 3.636509e-12
# 25 4.532433e-15 6.592625e-18 1.790600e-17 5.865462e-14
# 37 8.225917e-12 8.909007e-19 5.082756e-18 7.897758e-16
# 2  3.909235e-09 3.630384e-12 3.730264e-12 9.999995e-01

# look at the allele probabilities for the first five mice at marker "UNCHS000005"
apr$`1`[1:5,,"UNCHS000005"]
# A            B            C            D            E            F            G            H
# 1  4.699972e-09 4.681722e-09 5.145723e-04 4.994853e-01 4.699972e-09 4.999999e-01 5.328434e-08 1.210902e-07
# 13 9.906421e-10 9.793325e-10 6.046399e-08 6.046399e-08 9.906421e-10 4.999999e-01 5.320774e-08 4.999999e-01
# 25 9.365670e-09 9.329214e-09 5.000005e-01 4.999993e-01 9.365670e-09 1.361158e-11 3.114243e-11 2.419059e-07
# 37 7.130104e-04 9.985739e-01 1.396857e-08 1.396857e-08 7.130104e-04 7.716670e-11 4.402575e-10 5.768414e-09
# 2  1.954620e-09 1.932015e-09 1.208126e-07 1.208126e-07 1.954620e-09 1.815193e-12 1.865133e-12 9.999998e-01

########## perms.Rdata ############
# A permutation test establishes the statistical significance of a genome scan
# perms contains those values for n permutations (n=1000 in our case) for each phenotype
perms[1:5,] # look at first 5 perms 
#         fc       ob     ppi6    ppi6F
# [1,] 5.361236 4.874886 5.124117 5.075574
# [2,] 5.962797 6.418623 8.357625 8.254631
# [3,] 5.720598 6.397172 6.526800 7.003677
# [4,] 6.649054 5.802247 6.299134 6.448481
# [5,] 5.952524 5.665016 7.732061 7.249156

#extract values for specific thresholds for each phenotype
summary(perms, alpha=c(0.1, 0.05))
# LOD thresholds (1000 permutations)
#       fc   ob ppi6 ppi6F
# 0.1  7.22 7.23 7.22  7.22
# 0.05 7.70 7.64 7.67  7.66


#####################################################
          # Processing Scan Results 
#####################################################

####### Create a LOD plot for a single phenotype #############
# create a dataframe with .10 and .05 thesholds for the phenotypes
thr = data.frame(summary(perms,alpha=c(.10,.05)))
# thr       fc       ob     ppi6    ppi6F
# 0.1  7.217574 7.232223 7.220365 7.224730
# 0.05 7.701073 7.640308 7.669893 7.660169


# lets plot the fourth phenotype, ppi6F
#set up our plot area
par(mar=c(5.1, 4.1, 1.1, 1.1))

# to save the file to a png rather than view in the plot window in RStudio, uncomment the next line
#png(file=paste0(colnames(out)[4],".all.chr.",cM.or.Mbp,".png"), width = 1050, height=600)}

# find max y value for plotting
ymx = max(thr[,"ppi6F"])

# define map as genetic (cM) or physical (Mbp)
# we are using genetic map here
map=do$gmap

#plot the LOD scores all chromosomes
plot(out, map, lodcolumn="ppi6F", col="purple", ylim=c(0, ymx*1.02))

# add a title that has the phenotype, map type, and number of # subjects
title(main=paste0("ppi6F", "  ", sub="cM", "  ", length(which(!is.na(do$pheno[,"ppi6F"]))), " subjects"))

# add lines for the .05 and .10 thresholds 
abline(h=thr[2,"ppi6F"],col="red")
abline(h=thr[1,"ppi6F"], col="slategrey")

# add legend
legend("bottomleft", lwd=2, col="purple", "ppi6F", bg="gray90")

# if saving to png rather than outputting to plot window, uncomment next line.
# dev.off()


# since we have 4 phenotypes, rather than copy/paste the code above for each pheno,
# use a convenience function that will plot one plot for each phenotype
# The first will plot use the genetic (cM) map
plotLOD.all.chr(out, do, "cM", png=FALSE)
#The second will use the physical (Mbp) map
plotLOD.all.chr(out, do, "Mbp", png=FALSE)

#To save the above plots to a .png use png=TRUE in the function calls above

# We have a significant peak on chr13 for phenotype ppi6F
# Plot LOD plot for the peak for ppi6F on chromosome 13 only

# get max marker for that phenotype, map can be genetic or physical
# we're using genetic here
(ppi6F.chr13.max.marker.cM = max(out, lodcolumn="ppi6F", map=do$gmap))
# marker      chr      pos    ppi6F
# UNCHS036778  13 42.30464 8.220053

# get marker name
(ppi6F.marker=rownames(ppi6F.chr13.max.marker.cM))
#[1] "UNCHS036778"

# get allele probabilities for that marker
aprobs=getAlleleProbs(marker=ppi6F.marker,chr=13,apr=apr)
# get pheno data
pheno.dat=getPhenoData(do=do,pheno="ppi6F")
#write data to file
createProbsPhenoTable(probs=aprobs,pheno=pheno.dat,file="ppi6F.chr13.AlleleProbs.txt")

# for all phenotypes in out.Rdata, find all peaks above a given LOD threshold
# print the pheno, chr, low end of 95% CI, peak, high end of 95% CI and LOD
peaks.w.cis.Mbp = findPeaks(out=out,map=do$gmap,lod.thr=6)
# phenotype chromosome    ci_low position  ci_high      LOD
# 1         fc          6 39.797204 48.25264 48.77996 6.607959
# 11        ob          1 57.544957 57.87927 72.34658 6.334309
# 12      ppi6         13 41.614488 42.30464 42.86974 7.025436
# 13      ppi6         14  2.020275 20.87267 23.03369 6.128996
# 14     ppi6F         13 41.664778 42.30464 42.67592 8.220053

# write all the peaks with CIs out to a table with both locations in both Mbp and cM
write.table(peaks.w.cis.Mbp, file="peaks.with.CIs.Mbp.cM.txt", sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE, append=FALSE)
# do the same with units in cM 
peaks.w.cis.cM = findPeaks(out,do$gmap,lod.thr=5)
write.table(peaks.w.cis.cM, file="peaks.with.CIs.Mbp.cM.txt", sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE, append = TRUE)


# get the info on the 95% CI for ppi6F, chr13
(interval = data.frame(peaks.w.cis.Mbp[which(peaks.w.cis.Mbp$phenotype=="ppi6F" & peaks.w.cis.Mbp$chromosome==13),c("chromosome","ci_low","ci_high")]))
# chromosome   ci_low  ci_high
# 14         13 41.66478 42.67592


# write out table of genes in 95% CI
# requires DOQTL package
# install and load DOQTL
source("https://bioconductor.org/biocLite.R")
biocLite("DOQTL")
library(DOQTL)
(genes= getMGIGenesInCI("ppi6F",interval[1,"chromosome"],interval[1,"ci_low"],interval[1,"ci_high"],"gene","MGI"))

## phenotype x genotype plot
# infer genotype at a specific QTL position
# we'll use our chr 13 peak for ppi6F, locaated @ 42.30464 cM
g <- maxmarg(pr, do$gmap, chr=13, pos=42.30464, return_char=TRUE)
# plot genotype vs phenotype
# uncomment out 1st and 3rd lines below to save as .png
#png(file="genoxpheno.ppi6F.png")
plot_pxg(g, do$pheno[,"ppi6F"])
#dev.off()
```

Next Step [R/qtl2 Tutorial](https://github.com/Sethupathy-Lab/R-qtl2-pipeline/blob/master/Rqtl2Tutorial.Rmd)




 




