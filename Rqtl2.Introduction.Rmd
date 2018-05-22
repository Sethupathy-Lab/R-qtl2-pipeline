# R/qtl2-pipeline for Diversity Outbred(DO) mice
One of current protocols for QTL analysis of DO mice uses R/qtl2 written by [Karl Broman](http://kbroman.org/pages/about.html)

Reading the [documention ](http://kbroman.org/qtl2/docs.html)for R/qtl2 should be your starting point. 

There is also a good description of QTL mapping and some of the data structures produced in the qtl analysis in this [tutorial by Susan McClatchy ](https://smcclatchy.github.io/mapping/aio/). Working through this tutorial would give you a nice understanding of QTL mapping and the data structures in qtl2.

####Please note that, throughout this documentation, Rqtl2 and qtl2 are used interchangeably to refer to the R/qtl2 software package.  

There are 4 major parts to processing the data:  
1) Acquiring/Setting up the software  
2) Preparing the DO mouse Data  
3) Running the QTL scans  
4) Processing the scan output from the previous step  
All of these will be covered in separate documents.  


For those who would like a deeper understanding of R/qtl and QTL analysis in general, Karl has also written a book, "A Guide to QTL Mapping with R/qtl".  Although this book uses a previous version of the software (R/qtl vs R/qtl2) the QTL theory is the same and R/qtl2 is simply a reimplemenation of R/qtl.  A pdf copy the book is included in the R/qtl2 Tutorial which can found in the Sethupathy Lab directory on the CU server (/home/pr46_0001/cornell_tutorials/Rqtl2_Tutorial), which may be downloaded and read, but not distributed. The book can be purchased from [Amazon ](https://www.amazon.com/Guide-Mapping-Statistics-Biology-Health/dp/1461417082), [Springer ](https://www.springer.com/us/book/9780387921242?gclid=Cj0KCQjwre_XBRDVARIsAPf7zZitES-S1XgYyDUJOBDirplYflNJIwCtHPeQ_Twy0ZssOoLcI18CFbsaAl75EALw_wcB), and other online retailers.

From the introduction on the [Springer web page](https://www.springer.com/us/book/9780387921242?gclid=Cj0KCQjwre_XBRDVARIsAPf7zZitES-S1XgYyDUJOBDirplYflNJIwCtHPeQ_Twy0ZssOoLcI18CFbsaAl75EALw_wcB)...

"Quantitative trait locus (QTL) mapping is used to discover the genetic and molecular architecture underlying complex quantitative traits. It has important applications in agricultural, evolutionary, and biomedical research. R/qtl is an extensible, interactive environment for QTL mapping in experimental crosses. It is implemented as a package for the widely used open source statistical software R and contains a diverse array of QTL mapping methods, diagnostic tools for ensuring high-quality data, and facilities for the fit and exploration of multiple-QTL models, including QTL x QTL and QTL x environment interactions. This book is a comprehensive guide to the practice of QTL mapping and the use of R/qtl, including study design, data import and simulation, data diagnostics, interval mapping and generalizations, two-dimensional genome scans, and the consideration of complex multiple-QTL models. Two moderately challenging case studies illustrate QTL analysis in its entirety.

The book alternates between QTL mapping theory and examples illustrating the use of R/qtl. Novice readers will find detailed explanations of the important statistical concepts and, through the extensive software illustrations, will be able to apply these concepts in their own research. Experienced readers will find details on the underlying algorithms and the implementation of extensions to R/qtl. There are 150 figures, including 90 in full color.

Karl W. Broman is Professor in the Department of Biostatistics and Medical Informatics at the University of Wisconsin-Madison, and is the chief developer of R/qtl. Saunak Sen is Associate Professor in Residence in the Department of Epidemiology and Biostatistics and the Center for Bioinformatics and Molecular Biostatistics at the University of California, San Francisco."

Next Step
[Installing R/qtl2 software](https://github.com/Sethupathy-Lab/R-qtl2-pipeline/blob/master/Rqtl2.InstallingRqtl2.Rmd) 







