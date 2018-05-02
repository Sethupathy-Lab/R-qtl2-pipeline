# R/qtl2-pipeline for Diversity Outbred(DO) mice
One of current protocols for QTL analysis of DO mice uses R/qtl2 written by Karl Broman. His bio and contact information can be found here. http://kbroman.org/pages/about.html

Documentation for R/qtl2 can be found here and reading this documention should be your starting point.
http://kbroman.org/qtl2/docs.html

# Preparing DO mouse data
DO mice are a bit of special case and there are several steps for preparing the DO mouse data for R/qtl2.
A set of instructions for preparing the DO data can be found here.
http://kbroman.org/qtl2/pages/prep_do_data.html

## There are 4 basic steps to prepare the DO data
## 1) Download the Founder genotypes and SNP maps
With DO data, an important component of the R/qtl2 input files is the SNP genotypes for the eight founder lines. We also need the genetic and physical positions of the SNP markers.

Dan Gatti (https://www.jax.org/research-and-faculty/faculty/research-scientists/daniel-gatti) has prepared such data for the MegaMUGA and GigaMUGA arrays and made them available at ftp.jax.org/MUGA. The files with names starting MM_ are for the MegaMUGA array; the files with names starting GM_ are for the GigaMUGA array.

We’ve also placed these files at figshare: https://doi.org/10.6084/m9.figshare.c.3879694.

MM_primary_files.zip (https://doi.org/10.6084/m9.figshare.5404717.v1), containing all of the files for the MegaMUGA array
GM_primary_files.zip (https://doi.org/10.6084/m9.figshare.5404675.v1), containing all of the files for the GigaMUGA array
The files at ftp://ftp.jax.org/MUGA and in these zip files are all in .Rdata format, for R. For the GigaMUGA data, the key files include GM_snps.Rdata, containing SNP locations, and GM_geno.Rdata, containing SNP genotypes for the Collaborative Cross (CC) founder strains (which are also the founders of the DO), as single-letter codes N/C/T/A/G/H.

To use these data with R/qtl2, we need to reformat the files and recode the SNP genotypes (for example, as A/H/B/-). We’ve placed an R script to do this at the figshare site, parse_muga.R (https://doi.org/10.6084/m9.figshare.5405260). This uses the R/qtl2convert package, in particular find_consensus_geno() to determine consensus genotypes for each founder strain at each SNP from the multiple individuals that were genotyped, find_unique_geno() to determine the SNP alleles, and encode_geno() to encode the SNP genotypes as A/H/B/-.

The processed files are also at figshare:

MM_processed_files.zip (https://doi.org/10.6084/m9.figshare.5404750)
GM_processed_files.zip (https://doi.org/10.6084/m9.figshare.5404759)
MMnGM_processed_files.zip (https://doi.org/10.6084/m9.figshare.5404762)
The MM_ and GM_ files are for the MegaMUGA and GigaMUGA arrays, respectively. The MMnGM_ file (produced with the R script combine_MMnGM.R (https://doi.org/10.6084/m9.figshare.5405281.v1)) merges the SNPs for the two arrays, for DO projects that have some mice genotyped with the MegaMUGA array and others genotyped with the GigaMUGA array.

Each of these zip files contain a series of CSV files. For example, for the GigaMUGA data, we have:

GM_allelecodes.csv, the allele encodings for each SNP
GM_foundergeno*.csv, the CC founder genotypes, encoded as A/H/B/-, with one file for each chromosome
GM_gmap*.csv, the genetic map locations for the SNPS (in cM), with one file for each chromosome
GM_pmap*.csv, the physical map locations for the SNPS (in Mbp), with one file for each chromosome
GM_info.csv, information for all SNPs, including the genetic and physical locations and tier (an indication of SNP quality)
If your DO data includes just the GigaMUGA array, you’ll need the GM_processed_files.zip file. If you’ve used just the MegaMUGA array, you’ll need the MM_processed_files.zip file. If your project includes some MegaMUGA and some GigaMUGA arrays, you’ll need the MMnGM_processed_files.zip file.

Download the appropriate file and unzip it somewhere.









