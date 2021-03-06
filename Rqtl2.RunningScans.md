## 3) Running the QTL scans

Karl has a section in his user guide dedicated to [QTL scans in DO mice](http://kbroman.org/qtl2/assets/vignettes/user_guide.html#qtl_analysis_in_diversity_outbred_mice).  The script below is taken from that section of his guide and from discussions with Karl.  It would be prudent to read the entire guide, not just the DO portion, to better understand what each part of the scan is doing and how Rqtl2 works in general.

```
# reference your local R repository and load qtl2 library
.libPaths(c(.libPaths(), "/home/pr46_0001/shared/R_libs"))
library(qtl2)

# save your objects at each step, they will be 
# processed in later steps to produce plots etc.
# Also, if your process fails somewhere in the middle, you
# can just load the objects that are already saved to avoid 
# having to recreate them.  

# read in the DO data
do = read_cross2("qtl2.control.json")
save(do, file="do.Rdata")
#load("do.Rdata")

# calculate genotype probs
pr = calc_genoprob(do, error_prob = 0.002, cores = 0)
save(pr, file="pr.Rdata")

# convert genotype probs to allele probs
apr <- genoprob_to_alleleprob(pr)
save(apr, file="apr.Rdata")
#load("apr.Rdata")

# calculate kinship
k <- calc_kinship(apr, "loco")
save(k, file="k.Rdata")
#load("k.Rdata")

# create covar for sex
sex <- (do$covar$sex == "male")*1
names(sex) <- rownames(do$covar)
save(sex, file="sex.Rdata")
#load("sex.Rdata")

# genome scan with a linear mixed model (adjusting
# for a residual polygenic effect), with sex as an additive covariate.
out <- scan1(apr, do$pheno, k, sex, cores = 0)
save(out, file="out.Rdata")
#load("out.Rdata")

# permuations test
perms <- scan1perm(apr, do$pheno, k, n_perm=1000,  cores=0)
save(perms, file="perms.Rdata")
```
Next Step [Processing Scan Results](https://github.com/Sethupathy-Lab/R-qtl2-pipeline/blob/master/Rqtl2.ProcessingScanResults.md)
