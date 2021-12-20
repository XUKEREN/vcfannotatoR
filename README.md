# :dna: vcfannotatoR :dna: <!-- omit in toc -->
  [![HitCount](http://hits.dwyl.com/XUKEREN/vcfannotatoR.svg?style=flat-square)](http://hits.dwyl.com/XUKEREN/vcfannotatoR)
A variant annotation tool that parses vcf files and fetches variant information from the Ensembl Variant Effect Predictor (VEP) REST API.  

- [Introduction](#introduction)
- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Install vcfannotatoR from source](#install-vcfannotator-from-source)
- [How to use vcfannotatoR](#how-to-use-vcfannotator)
  - [Arguments](#arguments)
  - [Examples](#examples)
    - [Command line examples](#command-line-examples)
    - [Example inputs](#example-inputs)
    - [Example outputs](#example-outputs)
- [Contact](#contact)
## Introduction   

vcfannotatoR parses vcf files and annotates each variant in the vcf with the following information:
1. Type of variation (substitution, insertion, CNV, etc.).
2. Functional consequence (missense, silent, intergenic, etc.). If there are multiple effects, the variant will be annotated with the most deleterious consequence.
3. Sequence reading depth at each variant site.
4. Number of reads supporting the alternative allele.
5. Percentage of reads supporting the alternative allele versus those supporting the reference allele.
6. Allele frequency of variant (1000 genomes project) from Ensembl Variant Effect Predictor (VEP) REST API (API documentation is available here: http://grch37.rest.ensembl.org/documentation/info/vep_hgvs_post).
7. Additional annotations: rsid, minor_allele_freq, minor_allele, clinvar_significance, pubmed, transcript_id, gene_id, impact (a subjective classification of the severity of the variant consequence, based on agreement with SNPEff), gene_symbol (for each transcript), biotype, polyphen_prediction, and sift_prediction  

vcfannotatoR is open-source and available on [GitHub](https://github.com/XUKEREN/vcfannotatoR "GitHub vcfannotatoR page")

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
```

vcfannotatoR will check if the required R libraries are installed and automatically installs those not found. You can also use the following command in R to download the libraries before running vcfannotatoR:  


```r
# create a function to check dependencies and install packages if they are not already installed <!-- omit in toc -->
check_package <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages(package, dependencies = TRUE, repos='http://cran.us.r-project.org')
    library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
}

# check dependencies and install packages if they are not already installed <!-- omit in toc -->
packages <- c("optparse", "data.table", "tidyverse", "httr", "jsonlite", "xml2", "foreach", "doParallel")
load_packages <- lapply(packages, check_package)
```
**note:** "data.table" and "tidyverse" are used for data manipulation, "optparse" is used to catch the arguments from the command line, "foreach" and "doParallel" are used to perform for loops in parallel, and "jsonlite", "httr", and "xml2" are used to fetch variant consequences from the Ensembl VEP REST API.  

### Install vcfannotatoR from source

You can clone the GitHub repository:  
```bash
git clone https://github.com/XUKEREN/vcfannotatoR.git
cd vcfannotatoR
```
**vcfannotatoR.R** is ready for use in your command line.   
## How to use vcfannotatoR     
To run vcfannotatoR execute the vcfannotatoR.R script. This script catches the arguments from the command line and produces a tsv with variant annotations.   
### Arguments  

You can type the following in the command line to check available options:  
```bash
Rscript vcfannotatoR.R -h  
```  

Argument | Description | Required
------------ | ------------ | ------------
--input_dir | directory with the input vcf file | **Yes**
-I, --input | vcf file to be processed | **Yes**
-M, --getmeta | TRUE or FALSE, indicates whether to generate meta data tsv | No (default = FALSE)
-F, --getinfo | TRUE or FALSE, indicates whether to generate info fields tsv | No (default = FALSE)
-M, --getformat | TRUE or FALSE, indicates whether to generate format fields tsv | No (default = FALSE)
--totalread | format field total read depth | No (default = DP)
--refread | format field reference allele read depth | No (default = RO)
--altread | format field alternative allele read depth | No (default = AO)
--varianttype | format field variant type | No (default = TYPE)

### Examples   
#### Command line examples    
**Use case 1:** You can type the following in the command line, which will return a tsv file `Challenge_data.annotated.tsv` in the input directory:  
```bash
Rscript vcfannotatoR.R --input_dir ./data -I Challenge_data.vcf
```
An example output can be found here: [Challenge_data.annotated.tsv](/data/Challenge_data.annotated.tsv)   

**Use case 2:** You can type the following in the command line, which will return a complete set of output files in the input directory:     
```bash
Rscript vcfannotatoR.R --input_dir ./data -I Challenge_data.vcf --getmeta TRUE --getinfo TRUE --getformat TRUE 
```
**note:** User can define the format fields for total read depth, ref allele read depth, alt allele read depth, and variant type.  The default fields are DP, RO, AO, and TYPE. For VCF that does not have a variant type field, user should assign any available field to --varianttype. 

#### Example inputs   
A typical input vcf file is provided under `./data` [Challenge_data.vcf](data/Challenge_data.vcf)  
#### Example outputs  
A complete set of output files are provided under `./data`  
- [Challenge_data.annotated.tsv](/data/Challenge_data.annotated.tsv) has annotations for each variant     

        CHR	POS	REF	ALT	total_read_depth	TYPE	ref_read_depth	alt_read_depth	AD_alt_vs_ref	most_severe_consequence	variant_allele_freq.1kg	rsid	minor_allele_freq	minor_allele	clinvar_significance	pubmed	transcript_id	gene_id	impact	gene_symbol	biotype	polyphen_prediction	sift_prediction
        1	931393	G	T	4124	snp	4029	95	0.0235790518739141														
        1	935222	C	A	1134	snp	480	652	1.35833333333333	missense_variant	0.4938	rs2298214	0.4938	A	NA	32203549	ENST00000428771	ENSG00000188290	MODERATE	HES4	protein_coding	benign	tolerated_low_confidence
        1	1277533	T	C	786	snp	0	786	Inf	synonymous_variant	0.998	rs307362	0.002	T	NA	NA	ENST00000378888,ENST00000378891	ENSG00000107404,ENSG00000107404	LOW,LOW	DVL1,DVL1	protein_coding,protein_coding	NA,NA	NA,NA

- [Challenge_data.info_meta.tsv](/data/Challenge_data.info_meta.tsv) has meta data for the info fields  

        INFO	description
        NS	Number of samples with data
        DP	Total read depth at the locus
        DPB	Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype 

- [Challenge_data.format_meta.tsv](/data/Challenge_data.format_meta.tsv) has meta data for the format fields  

        FORMAT	description
        GT	Genotype
        GQ	Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype
        GL	Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy

- [Challenge_data.info.tsv](/data/Challenge_data.info.tsv) has all the info fields  

        CHROM	POS	ID	REF	ALT	QUAL	FILTER	AB	ABP	AC	AF	AN	AO	CIGAR	DP	DPB	DPRA	EPP	EPPR	GTI	LEN	MEANALT	MQM	MQMR	NS	NUMALT	ODDS	PAIRED	PAIREDR	PAO	PQA	PQR	PRO	QA	QR	RO	RPL	RPP	RPPR	RPR	RUN	SAF	SAP	SAR	SRF	SRP	SRR	TYPE
        1	931393	.	G	T	2.17938e-13	.	0	0	0	0	6	95	1X	4124	4124	0.999031	9.61615	316.776	0	1	1	59.7579	65.2274	2	1	591.29	0.989474	0.966741	0	0	0	0	3774	160284	4029	51	4.13032	101.278	44	1	40	8.15326	55	1663	269.369	2366	snp
        1	935222	.	C	A	16866.7	.	0.574956	58.3503	4	0.666667	6	652	1X	1134	1134	0	44.7878	169.779	0	1	2	70	70	2	1	53.367	0.990798	0.975	0	0	0	0	24492	19222	480	398	72.0711	214.077	254	1	28	1186.05	624	16	910.975	464	snp
        1	1277533	.	T	C	28168.6	.	0	0	6	1	6	786	1X	786	786	0	94.5216	0	0	1	1	70	0	2	1	307.075	0.977099	0	0	0	0	0	31532	0	0	474	75.5143	0	312	1	376	6.20397	410	0	0	0	snp

- [Challenge_data.format.tsv](/data/Challenge_data.format.tsv) has all the format fields with different samples in separate rows     

        CHROM	POS	ID	REF	ALT	Indiv	GT	GQ	DP	DPR	RO	QR	AO	QA
        1	931393	.	G	T	normal	0/0/0	132.995	2063	2063,0	2063	82063	0	0
        1	931393	.	G	T	vaf5	0/0/0	132.995	2061	2061,95	1966	78221	95	3774
        1	935222	.	C	A	normal	0/1/1	160.002	567	567,326	240	9611	326	12246

## Contact
Keren Xu (kerenxu@usc.edu)

https://xukeren.rbind.io/ 
