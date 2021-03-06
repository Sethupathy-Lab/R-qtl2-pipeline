## Running the R/qtl2 Tutorial
There is a complete tutorial on the CU servers in the Sethupathy Lab directory.  
It is called Rqt2_tutorial and can be found here...  
/home/pr46_0001/cornell_tutorials  
All the files needed to run it are contained in the Rqtl2_tutorial folder.  
To run the tutorial...  
#### make a reservation on a workstation  
#### login to your login node and then
```
# ssh to the workstation and cd to the /workdir
ssh your_user_name@workstation_name.tc.cornell.edu
# cd to the workdir on the workstation
cd /workdir

# copy the tutorial folder to /workdir
cp -R /home/pr46_0001/cornell_tutorials/Rqtl2_tutorial/ ./
```
At this point you should explore the contents of the tutorial folder.  There are several .R files that you will run to prepare the data.  Before you run them, use cat <filename> or nano <filename> to look at the content of the various .R files to see what each is doing.
To actually run the files you have two choices. The recommended approach for now is to run each one separately as described below.  Once you understand the process you can just use "sh PrepareDOData.sh" or "nohup sh PrepareDOData.sh &" to run them all.

```
# run the file to download the GigaMUGA Founder genotypes and SNP maps 
# interactively
Rscript downloadGMFiles.R
# in the background
nohup Rscript downloadGMFiles.R &

<!-- If you ran it interactively you should see something like...   -->

<!-- trying URL 'https://ndownloader.figshare.com/files/9311209' -->
<!-- Content type 'binary/octet-stream' length 5020437 bytes (4.8 MB) -->
<!-- ================================================== -->
<!-- downloaded 4.8 MB -->

<!-- [1] "checking to see if GM files are there" -->
<!-- total 16340 -->
<!-- drwxrwxr-x 2 wp244 wp244    4096 May 14 12:35 . -->
<!-- drwxrwxr-x 5 wp244 wp244    4096 May 14 12:35 .. -->
<!-- -rw-rw-r-- 1 wp244 wp244 2088395 May 14 12:35 GM_allelecodes.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  154210 May 14 12:35 GM_foundergeno10.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  180298 May 14 12:35 GM_foundergeno11.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  146671 May 14 12:35 GM_foundergeno12.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  149682 May 14 12:35 GM_foundergeno13.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  143293 May 14 12:35 GM_foundergeno14.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  129204 May 14 12:35 GM_foundergeno15.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  124118 May 14 12:35 GM_foundergeno16.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  123534 May 14 12:35 GM_foundergeno17.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  113222 May 14 12:35 GM_foundergeno18.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244   87933 May 14 12:35 GM_foundergeno19.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  237431 May 14 12:35 GM_foundergeno1.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  242278 May 14 12:35 GM_foundergeno2.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  178651 May 14 12:35 GM_foundergeno3.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  184695 May 14 12:35 GM_foundergeno4.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  184240 May 14 12:35 GM_foundergeno5.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  183276 May 14 12:35 GM_foundergeno6.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  178729 May 14 12:35 GM_foundergeno7.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  161263 May 14 12:35 GM_foundergeno8.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  166781 May 14 12:35 GM_foundergeno9.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244     677 May 14 12:35 GM_foundergenoM.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  113501 May 14 12:35 GM_foundergenoX.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244    1016 May 14 12:35 GM_foundergenoY.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  175880 May 14 12:35 GM_gmap10.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  205447 May 14 12:35 GM_gmap11.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  167256 May 14 12:35 GM_gmap12.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  170535 May 14 12:35 GM_gmap13.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  164065 May 14 12:35 GM_gmap14.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  147305 May 14 12:35 GM_gmap15.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  141458 May 14 12:35 GM_gmap16.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  140782 May 14 12:35 GM_gmap17.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  129123 May 14 12:35 GM_gmap18.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  100212 May 14 12:35 GM_gmap19.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  262548 May 14 12:35 GM_gmap1.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  267779 May 14 12:35 GM_gmap2.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  197576 May 14 12:35 GM_gmap3.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  204118 May 14 12:35 GM_gmap4.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  203618 May 14 12:35 GM_gmap5.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  202301 May 14 12:35 GM_gmap6.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  197280 May 14 12:35 GM_gmap7.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  178109 May 14 12:35 GM_gmap8.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  184135 May 14 12:35 GM_gmap9.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244     380 May 14 12:35 GM_gmapM.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  125290 May 14 12:35 GM_gmapX.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244     631 May 14 12:35 GM_gmapY.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244 4951855 May 14 12:35 GM_info.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  138350 May 14 12:35 GM_pmap10.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  161357 May 14 12:35 GM_pmap11.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  131132 May 14 12:35 GM_pmap12.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  133963 May 14 12:35 GM_pmap13.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  128262 May 14 12:35 GM_pmap14.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  114773 May 14 12:35 GM_pmap15.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  110082 May 14 12:35 GM_pmap16.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  109499 May 14 12:35 GM_pmap17.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  100438 May 14 12:35 GM_pmap18.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244   77933 May 14 12:35 GM_pmap19.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  206235 May 14 12:35 GM_pmap1.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  210165 May 14 12:35 GM_pmap2.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  154703 May 14 12:35 GM_pmap3.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  159989 May 14 12:35 GM_pmap4.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  159393 May 14 12:35 GM_pmap5.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  158528 May 14 12:35 GM_pmap6.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  155183 May 14 12:35 GM_pmap7.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  139036 May 14 12:35 GM_pmap8.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244  143644 May 14 12:35 GM_pmap9.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244     526 May 14 12:35 GM_pmapM.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244   98321 May 14 12:35 GM_pmapX.csv -->
<!-- -rw-rw-r-- 1 wp244 wp244     833 May 14 12:35 GM_pmapY.csv -->
<!-- [1] 0 -->
<!-- [1] "removingre the original zip file"   -->
```
If you ran it in the background your output will be in nohup.out
You can cat nohup.out periodically and see the output.  When the job is finished, the last line in nohup.out will be
[1]+  Done                    nohup Rscript downloadGMFiles.R
```
# look at contents of nohup.out
cat nohup.out

# run geneseek2qtl2.R interactively
Rscript geneseek2qtl2.R
# OR run it in the background
nohup Rscript geneseek2qtl2.R &

<!-- If you ran it interactively you should see output on the screen like.. -->

<!-- [your_user_name@workstation_name Rqtl2_tutorial]$ Rscript geneseek2qtl2.R -->
<!--  -File: GeneseekResults/Example_FinalReport.txt -->
<!--  -Reading data -->
<!-- Read 26073138 rows and 11 (of 11) columns from 1.147 GB file in 00:00:27 -->
<!--  -Reorganizing data -->
<!--  --Sample 10 of 182 -->
<!--  --Sample 20 of 182 -->
<!--  --Sample 30 of 182 -->
<!--  --Sample 40 of 182 -->
<!--  --Sample 50 of 182 -->
<!--  --Sample 60 of 182 -->
<!--  --Sample 70 of 182 -->
<!--  --Sample 80 of 182 -->
<!--  --Sample 90 of 182 -->
<!--  --Sample 100 of 182 -->
<!--  --Sample 110 of 182 -->
<!--  --Sample 120 of 182 -->
<!--  --Sample 130 of 182 -->
<!--  --Sample 140 of 182 -->
<!--  --Sample 150 of 182 -->
<!--  --Sample 160 of 182 -->
<!--  --Sample 170 of 182 -->
<!--  --Sample 180 of 182 -->
<!--  -Encode genotypes -->
<!--  -Grab X and Y intensities -->
<!--  --Sample 10 of 182 -->
<!--  --Sample 20 of 182 -->
<!--  --Sample 30 of 182 -->
<!--  --Sample 40 of 182 -->
<!--  --Sample 50 of 182 -->
<!--  --Sample 60 of 182 -->
<!--  --Sample 70 of 182 -->
<!--  --Sample 80 of 182 -->
<!--  --Sample 90 of 182 -->
<!--  --Sample 100 of 182 -->
<!--  --Sample 110 of 182 -->
<!--  --Sample 120 of 182 -->
<!--  --Sample 130 of 182 -->
<!--  --Sample 140 of 182 -->
<!--  --Sample 150 of 182 -->
<!--  --Sample 160 of 182 -->
<!--  --Sample 170 of 182 -->
<!--  --Sample 180 of 182 -->
<!--  -Writing X and Y intensities -->
<!--  -Writing genotypes -->



# List files in the directory to see the output of geneseek2qtl2.R
ls genos

<!-- there will be a bunch of csv files in your directory, e.g -->
<!-- example_geno10.csv    example_geno15.csv  example_geno1.csv  example_geno6.csv  example_genoX.csv -->
<!-- example_chrXint1.csv  example_geno11.csv  example_geno16.csv  example_geno2.csv  example_geno7.csv  example_genoY.csv -->
<!-- example_chrXint2.csv  example_geno12.csv  example_geno17.csv  example_geno3.csv  example_geno8.csv  geneseek2qtl2.R -->
<!-- example_chrYint1.csv  example_geno13.csv  example_geno18.csv  example_geno4.csv  example_geno9.csv example_chrYint2.csv      example_geno14.csv    example_geno19.csv  example_geno5.csv  example_genoM.csv -->
```
#### 3) Prepare the phenotype and covariate data as described in step 3 above.
There are example phenotype and covar files in the phenos folder  
You can look at them using ...
```
head phenos/example_phenos.csv
head phenos/example_covar.csv
```
#### 4) Prepare the control files by running the writeControlFiles.R script
```
# run writeControlFiles.R interactively
Rscript writeControlFiles.R
# OR run it in the background
nohup Rscript writeControlFiles.R &

<!-- This is the only output from this file. -->
<!-- [wp244@cbsulm15 Rqtl2_tutorial]$ Rscript writeControlFiles.R -->
<!-- Loading qtl2geno -->
<!-- Loading qtl2scan -->
<!-- Loading qtl2plot -->
<!-- If you don't see the 3 loading statements, then R is not recognizing your R_libs. -->
<!-- Check your code to make sure the statement .libPaths(c(.libPaths(),"/home/pr46_0001/shared/R_libs")) is the first line in the writeControlFiles.R file -->
```

#### Run the QTL scans
Running the scans takes a long time, long being a function of the how many samples in your data and the workstation you are on.  The medium memory machines may take several days. The largest large memory machines may crank it out in few hours.  Always reserve more time than you think you'll need.  Running the scans interactively is NOT recommended.  If your connection times out for some reason you will have to start over. 

```
# run the scans in the background
nohup Rscript runQTLScans.R &
```
You can periodically check your nohup.out file to see where the process is.  
Also be forwarned that you will receive **NO** error messages, so if your job dies, you will not know it unless you keep an eye on things.  

#### Process Scan Results
```
# list the .Rdata files that were output from running the scans
ls -l *.Rdata
<!-- -rw-rw-r-- 1 wp244 wp244  935677085 May 16 08:38 apr.Rdata -->
<!-- -rw-rw-r-- 1 wp244 wp244    9245922 May 16 08:38 do.Rdata -->
<!-- -rw-rw-r-- 1 wp244 wp244    4583749 May 16 08:38 k.Rdata -->
<!-- -rw-rw-r-- 1 wp244 wp244    2872531 May 16 08:38 out.Rdata -->
<!-- -rw-rw-r-- 1 wp244 wp244      22544 May 16 08:38 perms.Rdata -->
<!-- -rw-rw-r-- 1 wp244 wp244 4240158656 May 16 08:38 pr.Rdata -->
<!-- -rw-rw-r-- 1 wp244 wp244        443 May 16 08:38 sex.Rdata -->
```
Create a directory on your local computer and transfer the above .Rdata files from the scan results to that directory using an ftp program like Filezilla or some other file transfer utility like rsync.
Transfer the two .R files, processScanResults.R and processQtl2Functions.R to your local computer into the same directory as the .Rdata files. 
Finally, transfer the example_results file into this folder also. It contains some output from the code in processScanResults.R so that you can see what running the code should produce.  For best viewing, open the .txt file in Excel.

Open the processScanResults.R file in RStudio and change the setwd() command to point to the folder that you just created.  Run the code in a line by line fashion to learn about the data structures you just downloaded and also how to process the qtl2 results to get meaningful results.

There are several .png and .txt files that the above code puts out.  There are also examples of what these files look like in the example_results folder for comparison to your output. 

