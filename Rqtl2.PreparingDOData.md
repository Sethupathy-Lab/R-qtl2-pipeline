## 2) Preparing DO mouse data
DO mice are a bit of special case and there are several steps for preparing the DO mouse data for R/qtl2.  Because it is essential to have the data in the correct format, this page will cover DO data prep in detail.  Before attempting to run R/qtl2 it is imperative that you read and understand the [instructions for preparing the DO data](http://kbroman.org/qtl2/pages/prep_do_data.html), much of which is reproduced below.

### There are 4 basic steps to prepare the DO data
#### 1) Download the Founder genotypes and SNP maps
With DO data, an important component of the R/qtl2 input files is the SNP genotypes for the eight founder lines. The genetic and physical positions of the SNP markers are also required.

[Dan Gatti](https://www.jax.org/research-and-faculty/faculty/research-scientists/daniel-gatti) has prepared such data for the MegaMUGA and GigaMUGA arrays and made them available at [ftp.jax.org/MUGA](ftp://ftp.jax.org/MUGA). The files with names starting MM_ are for the MegaMUGA array; the files with names starting GM_ are for the GigaMUGA array.

These files are also available on [figshare ](https://doi.org/10.6084/m9.figshare.c.3879694).
[MM_primary_files.zip](https://doi.org/10.6084/m9.figshare.5404717.v1), contains all of the files for the MegaMUGA array.
[GM_primary_files.zip](https://doi.org/10.6084/m9.figshare.5404675.v1), contains all of the files for the GigaMUGA array.
The files at ftp.jax.org/MUGA and in the zip files on figshare are all in .Rdata format, for R. For the GigaMUGA data, the key files include GM_snps.Rdata, containing SNP locations, and GM_geno.Rdata, containing SNP genotypes for the Collaborative Cross (CC) founder strains (which are also the founders of the DO), as single-letter codes N/C/T/A/G/H.

To use the JAX and figshare files above with R/qtl2, the files must be reformated and recoded for the SNP genotypes (for example, as A/H/B/-). There is an R script to do this at the figshare site, [parse_muga.R](https://doi.org/10.6084/m9.figshare.5405260). This uses the R/qtl2convert package, in particular find_consensus_geno(), to determine consensus genotypes for each founder strain at each SNP from the multiple individuals that were genotyped, find_unique_geno() to determine the SNP alleles, and encode_geno() to encode the SNP genotypes as A/H/B/-.

**To avoid having to reformat/recode files, there is also a set of pre-processed files at figshare: **  
[MM_processed_files.zip](https://doi.org/10.6084/m9.figshare.5404750)  
[GM_processed_files.zip](https://doi.org/10.6084/m9.figshare.5404759)  
[MMnGM_processed_files.zip](https://doi.org/10.6084/m9.figshare.5404762)  

The MM_ and GM_ files are for the MegaMUGA and GigaMUGA arrays, respectively. The MMnGM_ file (produced with the R script [combine_MMnGM.R](https://doi.org/10.6084/m9.figshare.5405281.v1)) merges the SNPs for the two arrays, for DO projects that have some mice genotyped with the MegaMUGA array and others genotyped with the GigaMUGA array.

Each of these zip files contain a series of CSV files. For example, for the GigaMUGA data, we have:
GM_allelecodes.csv, the allele encodings for each SNP; GM_foundergeno\*.csv, the CC founder genotypes, encoded as A/H/B/-, with one file for each chromosome; GM_gmap\*.csv, the genetic map locations for the SNPS (in cM), with one file for each chromosome; GM_pmap\*.csv, the physical map locations for the SNPS (in Mbp), with one file for each chromosome; GM_info.csv, information for all SNPs, including the genetic and physical locations and tier (an indication of SNP quality)

If your DO data includes just the GigaMUGA array, you’ll need the GM_processed_files.zip file. If you’ve used just the MegaMUGA array, you’ll need the MM_processed_files.zip file. If your project includes some MegaMUGA and some GigaMUGA arrays, you’ll need the MMnGM_processed_files.zip file.

**Download the appropriate file to current working dir and unzip it. **  
To avoid having to reformat/recode files, it's easiest to download a set of pre-processed files at figshare. Assuming we have genotyped using the GigaMUGA array, download [GM_processed_files.zip ](https://doi.org/10.6084/m9.figshare.5404759) and unzip the files into a directory called "GM". 

Get the figshare URL for the files by navigating to the link above, right or control click on the Download button and select Copy Link Address.  Enter the link address,
(https://ndownloader.figshare.com/files/9311209), in the code below at url <-, as shown.

The code below will download the files and leave you with a folder called GM in the current directory, with all the necessary files in it. 
This code assumes that you are logged onto a workstation and are in /workdir/Rqtl2_tutorial
```
# To run this interactively,
# if not already started, start R 
/programs/R-3.4.2/bin/R

# in R
# file path to Rqtl2_tutorial
setwd("/workdir/Rqtl2_tutorial")
# download GM files from figshare URL
url <- 'https://ndownloader.figshare.com/files/9311209'
# download to current dir with file ext = .zip
path1 <- tempfile(tmpdir=getwd(),fileext = ".zip")
download.file(url, path1, mode="wb")
#unzip it
unzip(path1)
print("check to see if GM files are there")"
print(system("ls -al GM"))
# optional...you can remove the original zip file now
system("rm *.zip")
```

#### 2) Encoding the DO genotypes
GeneSeek will provide a zip file that contains a series of files. **FinalReport.txt** file, which is the biggest of them, is the only file needed.  It contains one line for each SNP and for each sample. For 200 mice, it’ll be about 1 GB, and maybe a quarter of that size when compressed as a zip file.

Provided is an example R script, [geneseek2qtl2.R](http://kbroman.org/qtl2/assets/geneseek2qtl2.R), to convert the FinalReport.txt file into what’s needed for R/qtl2. It does two things: encodes the genotypes as A/H/B/- and writes them to a series of CSV files, one per chromosome, and extracts the SNP array intensities for the SNPs on the X and Y chromosome, which is useful for verifying the sex of the mice.
To use this script, you’ll need to edit three lines near the top:

**codefile** defines the path to the GM_allelecodes.csv file (or MM_ or MMnGM_ version, if you’re using MegaMUGA or both arrays in your project)

**ifiles** defines the path to your FinalReport.txt file. This can be a vector of such file paths, if your genotyping was performed in batches

**ostem** defines the path and file “stem” for the output files. For example, if you use ostem <- "Data/forqtl2", the output files will be placed in the Data subdirectory and will have names like forqtl2_geno1.csv.

An issue you may need to contend with is a possible recoding of the sample identifiers. For example, you may have one batch where the samples are labeled like DO-146 and another where they’re labeled simply 146 and another where they’re labeled AA-DO-146. Search for the following line, and do some reformatting of the sample IDs there.  
**#Note: may need to revise the IDs in the second column**  
You’ll need to have the [R/qtl2convert package](https://github.com/rqtl/qtl2convert) installed. And note that, since this script reads in the full data into memory, you’ll need a computer with appreciable RAM.

The result of the [geneseek2qtl2.R](http://kbroman.org/qtl2/assets/geneseek2qtl2.R) script will be a series of CSV files containing the re-coded genotypes, with one file per chromosome, plus two files containing the SNP array intensities for the X and Y chromosomes (one file per allele per chromosome).  
```
# This is what the geneseek2qtl2.R file should look like
# reference your R repository
#.libPaths(c(.libPaths(),"~/R_libs"))
# or reference the shared R_libs in the Sethupathy Lab directory
.libPaths(c(.libPaths(),"/home/pr46_0001/shared/R_libs"))
# set you working directory to the Rqtl2_tutorial directory
setwd("/path/to your/directory")

# convert GeneSeek FinalReport files to format for R/qtl2
# requires the qtl2convert package, https://github.com/rqtl/qtl2convert
#
# - creates one genotype CSV file for each chromosome
#
# - also creates 4 files containing the two channels of SNP intensities for markers on the X and Y chr
#   (these are useful for verifying the sex of the mice)

# file containing allele codes for GigaMUGA data
#   - from GM_processed_files.zip, https://doi.org/10.6084/m9.figshare.5404759
codefile <- "GM/GM_allelecodes.csv"

# input files with GigaMUGA genotypes
#  - can be a single file or a vector of multiple files
#  - if samples appear in multiple files, the genotypes in later files
#    will be used in place of genotypes in earlier files
#  - files can be gzipped (".gz" extension)
ifiles <- "GeneseekResults/Example_FinalReport.txt"

# file "stem" for output files
# output files will be like "example_geno19.csv"
ostem <- "example"

##############################
# define a couple of functions
##############################
# simple version of data.table::fread()
myfread <- function(filename) data.table::fread(filename, data.table=FALSE)

# cbind, replacing matching columns with second set and adding unique ones
cbind_smother <-
  function(mat1, mat2)
  {
      cn1 <- colnames(mat1)
      cn2 <- colnames(mat2)
      m <- (cn2 %in% cn1)
      if(any(m)) {
          mat1[,cn2[m]] <- mat2[,cn2[m],drop=FALSE]
          if(any(!m)) {
              mat1 <- cbind(mat1, mat2[,cn2[!m],drop=FALSE])
          }
      }
      else {
          mat1 <- cbind(mat1, mat2)
      }

      mat1
  }
##############################

# read genotype codes
codes <- myfread(codefile)

full_geno <- NULL
cXint1 <- cXint2 <- NULL
cYint1 <- cYint2 <- NULL

for(ifile in ifiles) {
    cat(" -File:", ifile, "\n")
    rezip <- FALSE
    if(!file.exists(ifile)) {
        cat(" -Unzipping file\n")
        system(paste("gunzip", ifile))
        rezip <- TRUE
    }

    cat(" -Reading data\n")
    g <- myfread(ifile)
    # subset to the markers in the codes object
    g <- g[g[,1] %in% codes[,1],]

    # NOTE: may need to revise the IDs in the 2nd column
    # We're assuming the IDs are numbers and creating revised IDs as DO-###
    samples <- unique(g[,2])

    # matrix to contain the genotypes
    geno <- matrix(nrow=nrow(codes), ncol=length(samples))
    dimnames(geno) <- list(codes[,1], samples)

    # fill in matrix
    cat(" -Reorganizing data\n")
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (g[,2]==samples[i])
        geno[g[wh,1],i] <- paste0(g[wh,3], g[wh,4])
    }

    cat(" -Encode genotypes\n")
    geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))

    if(is.null(full_geno)) {
        full_geno <- geno
    } else {
        # if any columns in both, use those from second set
        full_geno <- cbind_smother(full_geno, geno)
    }

    # grab X and Y intensities
    cat(" -Grab X and Y intensities\n")
    gX <- g[g[,1] %in% codes[codes$chr=="X",1],]
    gY <- g[g[,1] %in% codes[codes$chr=="Y",1],]
    cX1 <- matrix(nrow=sum(codes$chr=="X"),
                  ncol=length(samples))
    dimnames(cX1) <- list(codes[codes$chr=="X",1], samples)
    cX2 <- cX1
    cY1 <- matrix(nrow=sum(codes$chr=="Y"),
                  ncol=length(samples))
    dimnames(cY1) <- list(codes[codes$chr=="Y",1], samples)
    cY2 <- cY1
    for(i in seq(along=samples)) {
        if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
        wh <- (gX[,2]==samples[i])
        cX1[gX[wh,1],i] <- gX$X[wh]
        cX2[gX[wh,1],i] <- gX$Y[wh]

        wh <- (gY[,2]==samples[i])
        cY1[gY[wh,1],i] <- gY$X[wh]
        cY2[gY[wh,1],i] <- gY$Y[wh]
    }
    if(is.null(cXint1)) {
        cXint1 <- cX1
        cXint2 <- cX2
        cYint1 <- cY1
        cYint2 <- cY2
    } else {
        # if any columns in both, use those from second set
        cXint1 <- cbind_smother(cXint1, cX1)
        cXint2 <- cbind_smother(cXint2, cX2)
        cYint1 <- cbind_smother(cYint1, cY1)
        cYint2 <- cbind_smother(cYint2, cY2)
    }

    if(rezip) {
        cat(" -Rezipping file\n")
        system(paste("gzip", ifile))
    }
}

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(cbind(marker=rownames(cXint1), cXint1),
                       paste0(ostem, "_chrXint1.csv"),
                       paste(ostem, "X chr intensities, channel 1"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cXint2), cXint2),
                       paste0(ostem, "_chrXint2.csv"),
                       paste(ostem, "X chr intensities, channel 2"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint1), cYint1),
                       paste0(ostem, "_chrYint1.csv"),
                       paste(ostem, "Y chr intensities, channel 1"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint2), cYint2),
                       paste0(ostem, "_chrYint2.csv"),
                       paste(ostem, "Y chr intensities, channel 2"),
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chr in c(1:19,"X","Y","M")) {
    mar <- codes[codes$chr==chr,1]
    g <- full_geno[mar,]
    qtl2convert::write2csv(cbind(marker=rownames(g), g),
                           paste0(ostem, "_geno", chr, ".csv"),
                           paste0(ostem, " genotypes for chr ", chr),
                           overwrite=TRUE)
}
```
#### 3) Preparing the phenotype and covariate data
You next need to prepare files with the phenotype data (with strictly numeric phenotypes) and covariate data (which can include non-numeric variables). Each should be arranged with samples as rows and variables as columns, and with the first column being the sample IDs, exactly as used in the genotype data.

For examples of these files, see the [qtl2data repository ](https://github.com/rqtl/qtl2data), in particular [do_pheno.csv ](https://github.com/rqtl/qtl2data/blob/master/DO_Gatti2014/do_pheno.csv) and [do_covar.csv ](https://github.com/rqtl/qtl2data/blob/master/DO_Gatti2014/do_covar.csv), for the [Gatti et al. (2014) ](http://www.g3journal.org/content/ggg/4/9/1623.full.pdf) data.

The covariate data needs to contain a sex column, plus a column giving the generation number (e.g., ngen) for each DO animal.

For help in determining the generation numbers, see the file [DO_generation_dates.csv ](http://kbroman.org/qtl2/assets/DO_generation_dates.csv) which lists the generation number for different dates of distribution of mice from the Jackson Lab.  

#### **Note that one assumption of qtl2 is that the phenotype data is normally distributed.  This is often not the case.  You should always examine the distributions of your phenotype data first and use an appropriate data transformation if it is not normally distributed.  

#### 4) Preparing the control file  
Finally, you need to prepare a control file which includes details of the file names and the various encodings and variable names.

Use the write_control_file() function as indicated below to make the control file. It’s a bit complicated, because there’s a lot of stuff to specify, and some of the pieces are sort of confusing.
```
# make sure you reference your local R_libs where qtl2 was installed
.libPaths(c(.libPaths(), "~/R_libs"))

library(qtl2)
chr <- c(1:19, "X")
write_control_file(# change forqtl2 to something meaningful for you 
                   "forqtl2.json", 
                   crosstype="do", #diversity outbred
                   description="My awesome DO project", # change to something meaningful for you
                   # the following commands with "GM/GM.."
                   # should point to the files you downlaoded
                   # in step 1
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("forqtl2_geno", chr, ".csv"), #forqtl should match **ostem** 
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="my_pheno.csv", # The forqtl2 portion of these two lines should 
                   covar_file="my_covar.csv", # match files created in Step 3, phenotype and covariate data
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen")  
```
The first three arguments are the name of the file to be created (with extension either .json or .yaml), the cross type ("do" for DO mice), and a description.

The argument founder_geno_file gives the name of the file that contains the founder genotypes, here a vector of file names since we have the genotypes split by chromosome. Then founder_geno_transposed=TRUE indicates that we have the founder strains as columns and the markers as rows, rather than the default orientation with strains as the rows.

The arguments gmap_file and pmap_file indicate the names of the files that contain the genetic (gmap) and physical (pmap) marker positions, respectively. These are again each a vector of file names, since the files are split by chromosome.

The argument geno_file indicates the name of the file with the DO genotype data (again, a vector of file names as they’re split by chromosome), and we use geno_tranposed=TRUE because the [geneseek2qtl2.R](http://kbroman.org/qtl2/assets/geneseek2qtl2.R) script discussed above writes the CSV files with individuals as columns and markers as rows.

The argument geno_codes indicates the encodings for both the DO genotypes and the founder genotypes. This may be a bit confusing, but we want to make a list with the names being the codes used in the data and the values being the integers 1, 2, and 3. So geno_codes=list(A=1, H=2, B=3).

The argument xchr indicates the chromosome identifier for the X chromosome.

The pheno_file and covar_file arguments indicate the names of the files with the phenotype and covariate data, respectively. We then provide some details about specific covariates that are important: sex (via sex_covar) plus sex_codes to explain the encodings, with the names F and M being the codes used in the data and the values "Female" and "Male" indicating which sexes the codes correspond to. And finally crossinfo_covar indicates the name of the covariate column that contains the cross information ("ngen" for number of generations, for cross type "do").

It’s sometimes useful to create chromosome-specific control files, for the case that you want to just load one chromosome worth of data. You can do this by first defining chr to be a single value, e.g. chr <- "1". You can then use exactly the same write_control_file() code, though perhaps changing the first line, with the name of the file to be produced. For example:
```
# make sure you reference your local R_libs where qtl2 was installed
.libPaths(c(.libPaths(), "~/R_libs"))

library(qtl2)
chr <- "1"
write_control_file(# change forqtl2 to something meaningful for you 
                   "forqtl2_chr1.json", 
                   crosstype="do",
                   description="My awesome DO project", # change to something meaningful for you
                   # the following commands with "GM/GM.."
                   # should point to the files you downlaoded
                   # in step 1
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("forqtl2_geno", chr, ".csv"), # forqtl2 should match **ostem**
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="my_pheno.csv", # these names should match pheno and covariate files created
                   covar_file="my_covar.csv", # in step 3 above
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen")  
```

Next Step [Running the Scans](https://github.com/Sethupathy-Lab/R-qtl2-pipeline/blob/master/Rqtl2.RunningScans.Rmd) 
