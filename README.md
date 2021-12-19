# vcfannotatoR
A variant annotation tool parses vcf files and fetches variant information from the Ensembl Variant Effect Predictor (VEP) REST API.  

- [vcfannotatoR](#vcfannotator)
  - [Introduction](#introduction)
  - [Installation](#installation)
    - [Dependencies](#dependencies)
- [check dependencies and install packages if they are not already installed](#check-dependencies-and-install-packages-if-they-are-not-already-installed)
    - [Install vcfannotatoR from source](#install-vcfannotator-from-source)
  - [How to use vcfannotatoR](#how-to-use-vcfannotator)
    - [Arguments](#arguments)
    - [Examples](#examples)
## Introduction   

vcfannotatoR parses vcf files and annotate each variant in the vcf with the following information:
1. Type of variation (substitution, insertion, CNV, etc.).
2. Functional consequence (missense, silent, intergenic, etc.). If there are multiple effects, the variant will be annotated with the most deleterious consequence.
3. Sequence reading depth at each variant site.
4. Number of reads supporting the alternative allele.
5. Percentage of reads supporting the alternative allele versus those supporting the reference allele.
6. Allele frequency of variant (1000 genomes project) from Ensembl Variant Effect Predictor (VEP) REST API (API documentation is available here: http://grch37.rest.ensembl.org/documentation/info/vep_hgvs_post).
7. Additional annotations: rsid, minor_allele_freq, minor_allele, clinvar significance, pubmed, transcript_id, gene_id, impact, gene_symbol, biotype, polyphen_prediction, and sift_prediction

## Installation

### Dependencies  

vcfannotatoR needs the following:
- **R** (tested on R version 4.1.1)
- **An internet connection**
- **The following R libraries:** (The number is the version tested during development)
```r 
    data.table (1.14.2)     tidyverse (1.3.1)       optparse (1.7.1)         
    foreach (1.5.1)         doParallel (1.0.16)     jsonlite (1.7.2)
    httr (1.4.2)            xml2 (1.3.3)            
````

vcfannotatoR will check if the required R libraries are installed and automatically installs those not found. You can also use the following command in R to download the libraries before running vcfannotatoR:  


```r
# create a function to check dependencies and install packages if they are not already installed <!-- omit in toc -->
check_package <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages(package, dependencies = TRUE, repos='http://cran.us.r-project.org')
    library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
}

# check dependencies and install packages if they are not already installed  
packages <- c("optparse", "data.table", "tidyverse", "httr", "jsonlite", "xml2", "foreach", "doParallel")
load_packages <- lapply(packages, check_package)
````
**note:** "data.table" and "tidyverse" are used for data manipulation, "optparse" is used to catch the arguments from the command line, "foreach" and "doParallel" are used to perform for loops in parallel, and "jsonlite", "httr", and "xml2" are used to fetch variant consequences from the Ensembl VEP REST API.  

### Install vcfannotatoR from source

You can clone the GitHub repository:  
```bash
git clone https://github.com/XUKEREN/vcfannotatoR.git
````
**vcfannotatoR.R** is ready for use in your command line.   
## How to use vcfannotatoR     
To run vcfannotatoR execute the vcfannotatoR.R script. This script catches the arguments from the command line and produces a tsv with variant annotations.   
### Arguments  

You can type the following in the command line to check options:  
```bash
Rscript vcfannotatoR.R -h  
````  

Argument | Description | Required
------------ | ------------ | ------------
--input_dir | directory with the input *vcf* file | **Yes**
-I, --input | vcf file to be processed | **Yes**
-M, --getmeta | TRUE or FALSE, indicates whether to generate meta data tsv | No (default = FALSE)
-F, --getinfo | TRUE or FALSE, indicates whether to generate info fields tsv | No (default = FALSE)
-M, --getformat | TRUE or FALSE, indicates whether to generate format fields tsv | No (default = FALSE)
### Examples   
You can type the following in the command line, which will return a tsv file `Challenge_data.annotated.tsv` in the input directory:  
```bash
Rscript vcfannotatoR.R --input_dir ./data -I Challenge_data.vcf
````
An example output can be found here: [Challenge_data.annotated.tsv]()   

You can type the following in the command line, which will return a complete set of output files in the input directory:     
```bash
Rscript vcfannotatoR.R --input_dir ./data -I Challenge_data.vcf --getmeta TRUE --getinfo TRUE --getformat TRUE 
````
Example outputs: 
- [Challenge_data.annotated.tsv](/data/Challenge_data.annotated.tsv) has annotations for each variant 
- [Challenge_data.info_meta.tsv](/data/Challenge_data.info_meta.tsv) has meta data for the info fields  
- [Challenge_data.format_meta.tsv](/data/Challenge_data.format_meta.tsv) has meta data for the format fields  
- [Challenge_data.info.tsv](/data/Challenge_data.info.tsv) has all the info fields  
- [Challenge_data.format.tsv](/data/Challenge_data.format.tsv) has all the format fields   