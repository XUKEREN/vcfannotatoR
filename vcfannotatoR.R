################################################################################
#
#   File name: vcfannotatoR.R
#
#   Author: Keren Xu ( kerenxu@usc.edu )
#   PhD Candidate in Epidemiology
#   Center for Genetic Epidemiology
#   University of Southern California
#   https://xukeren.rbind.io/
#
################################################################################

################################################################################
#
#   Description: A variant annotation tool. This script parses vcf files and annotates each variant in the vcf with the following information:
#   1. Type of variation (substitution, insertion, CNV, etc.).
#   2. Functional consequence (missense, silent, intergenic, etc.). If there are multiple effects, the variant will be annotated with the most deleterious consequence.
#   3. Sequence reading depth at each variant site.
#   4. Number of reads supporting the alternative allele.
#   5. Percentage of reads supporting the alternative allele versus those supporting the reference allele.
#   6. Allele frequency of variant (1000 genomes project) from Ensembl Variant Effect Predictor (VEP) REST API (API documentation is available here: http://grch37.rest.ensembl.org/documentation/info/vep_hgvs_post).
#   7. Additional annotations: rsid, minor_allele_freq, minor_allele, clinvar_significance, pubmed, transcript_id, gene_id, impact, gene_symbol, biotype, polyphen_prediction, sift_prediction
#   Command line use example: Rscript vcfannotatoR.R --input_dir ./data -I Challenge_data.vcf
#   Command line use example that returns all the output files: Rscript vcfannotatoR.R --input_dir ./data -I Challenge_data.vcf --getmeta TRUE --getinfo TRUE --getformat TRUE
#   
################################################################################

# clear workspace
rm(list = ls())

# ===============================================================================
#    Load libraries
# ===============================================================================

# create a function to check dependencies and install packages if they are not already installed
check_package <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    install.packages(package, dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
}

# check dependencies and install packages if they are not already installed
packages <- c("optparse", "data.table", "tidyverse", "httr", "jsonlite", "xml2", "foreach", "doParallel")
load_packages <- lapply(packages, check_package)

# ===============================================================================
#    Catch the arguments from the command line
# ===============================================================================

# create command line option lists
option_list <- list(
  make_option("--input_dir", type = "character", default = NA, 
              help = "directory with the input vcf"),
  make_option(c("-I", "--input"), type = "character", default = NA, 
              help = "vcf file to be processed", metavar = "character"),
  make_option(c("-M", "--getmeta"), type = "logical", default = FALSE, 
              help = "output meta data", metavar = "logical"),
  make_option(c("-F", "--getinfo"), type = "logical", default = FALSE, 
              help = "output info fields", metavar = "logical"),
  make_option(c("-G", "--getformat"), type = "logical", default = FALSE, 
              help = "output format fields", metavar = "logical"),
  make_option("--totalread", type = "character", default = "DP", 
              help = "format field total read depth", metavar = "character"),
  make_option("--refread", type = "character", default = "RO", 
              help = "format field reference allele read depth", metavar = "character"),
  make_option("--altread", type = "character", default = "AO", 
              help = "format field alternative allele read depth", metavar = "character"),
  make_option("--varianttype", type = "character", default = "TYPE", 
              help = "format field variant type", metavar = "character")
)

# parses command line options 
opt <- parse_args(OptionParser(option_list = option_list))

# ===============================================================================
#    Read input vcf
# ===============================================================================

# check if the vcf name or the vcf input dir are provided by the user
if (is.na(opt$input_dir) || is.na(opt$input)) {
  message("\n\nPlease provide the input vcf name and the input dir\n\n")
  q() # quit if the vcf name or the vcf input dir are not provided
}

# get the name, dir, basename of the input vcf and use the basename to generate the file name of the output tsv
input_vcf <- paste(opt$input_dir, opt$input, sep = "/")
basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", opt$input)
output_tsv <- paste(opt$input_dir, paste0(basename, ".annotated.tsv"), sep = "/")

# read input vcf to a data.table without keeping the metadata
my_vcf <- read.table(input_vcf, stringsAsFactors = F, header = F, sep = "\t")

# read input vcf to text lines keeping the metadata
strings <- readLines(input_vcf)

# ===============================================================================
#    Create functions
# ===============================================================================

# get colnames for the vcf `my_vcf`
get_header <- function(string) {
  pattern <- "^#CHROM"
  header_line_index <- min(grep(pattern, string))
  header_line <- string[header_line_index]
  my_header <- str_split(header_line, "\t") %>% unlist()
  return(my_header)
}

# convert metadata into data.table
get_metadata <- function(pattern, string) {
  meta_index <- grep(pattern, string)
  meta_line <- string[meta_index]
  df_meta <-
    meta_line %>%
    as_tibble(.name_repair = "minimal") %>%
    separate(value, c("drop", "keep"), sep = "=<ID=") %>%
    select(-"drop") %>%
    separate(keep, c(sub("\\^##", "", pattern), "keep"), sep = ",Number") %>%
    separate(keep, c("drop", "keep"), "Description=\"") %>%
    select(-"drop") %>%
    separate(keep, c("description", "drop"), "\">") %>%
    select(-"drop")
  return(df_meta)
}

# get info field (variant level information) 
## split the info field, e.g., remove string "DP=" to only the value and set DP as the column name
get_info_field_splitted <- function(x) {
  x %>%
    as_tibble(.name_repair = "minimal") %>%
    separate(value, c("colnames", "value"), sep = "=") %>%
    pivot_wider(names_from = colnames, values_from = value)
}
## create a data frame with each info field in a separate column
get_info <- function(vcf) {
  info_field_splitted <- vcf$INFO %>%
    str_split(., ";") %>%
    map(get_info_field_splitted) %>%
    rbindlist()
  variant_level_columns <- c(1:(which(colnames(vcf) == "INFO") - 1))
  df_vcf_info <- cbind(vcf[, variant_level_columns], info_field_splitted)
  return(df_vcf_info)
}

# a function that can extract a single info field  
get_single_info <- function(string, info) {
  pattern <- paste0("^", info, "=")
  single_info_index <- grep(pattern, string)
  single_info_line <- string[single_info_index]
  single_info_line %>% sub(paste0(info, "="), "", .)
}

# get format field (individual level information)
get_format <- function(vcf) {
  my_header_format <- str_split(vcf[1, which(colnames(vcf) == "FORMAT")], ":") %>% unlist()
  variant_level_columns_short <- c(1:which(colnames(vcf) == "ALT"))
  individual_level_columns <- c((which(colnames(vcf) == "FORMAT") + 1):dim(vcf)[2])
  df_vcf_format <-
    vcf[c(variant_level_columns_short, individual_level_columns)] %>%
    pivot_longer(
      cols = names(vcf[individual_level_columns]),
      names_to = "Indiv",
      values_to = "format_field"
    ) %>%
    separate(col = format_field, into = my_header_format, sep = ":")
  return(df_vcf_format)
}

# ===============================================================================
#    Annotate vcf
# ===============================================================================

# use function get_header to assign my_vcf with colnames
colnames(my_vcf) <- get_header(strings)

# use function get_single_info to extract four columns: reading depth, type of the variants, reference allele observation count and alternate allele observation count  
df_anno <- c(opt$totalread, opt$varianttype, opt$refread, opt$altread) %>%
  map(get_single_info, string = my_vcf$INFO %>% str_split(., ";") %>% unlist()) %>%
  do.call(cbind, .) %>%
  as_tibble(.name_repair = "minimal")
colnames(df_anno) <- c("DP", "TYPE", "RO", "AO")

# merge with variant level columns  
variant_level_columns <- c(1:(which(colnames(my_vcf) == "INFO") - 1))
df_anno <- cbind(my_vcf[, variant_level_columns], df_anno)

# split multi-allelic variant calls into separate lines
df_anno <- df_anno %>% separate_rows(., ALT, TYPE, AO, sep = ",", convert = TRUE)

# calculate percentage of reads supporting the alternative allele versus those supporting reference allele.
df_anno <- df_anno %>%
  mutate(AO = as.numeric(AO), RO = as.numeric(RO)) %>%
  mutate(AD_alt_vs_ref = AO / RO)

# create HGVS notations that will be fed into the Ensembl VEP REST API
hgvs_notations <- df_anno %>%
  mutate(hgvs_notation = paste0(`#CHROM`, ":g.", POS, REF, ">", ALT)) %>%
  pull(hgvs_notation)

# fetch variant consequences for multiple HGVS notations
server <- "http://grch37.rest.ensembl.org"
ext <- "/vep/human/hgvs"

# only can request 300 each time ... so set up iterations
iterations <- round(length(hgvs_notations) / 300) + 1

# use foreach to perform fetch request in a loop in parallel
# setup parallel backend to use many processors
cores <- detectCores()
# not to overload your computer
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)
message(paste0("\n\nfetching ", length(hgvs_notations), " variants... It might take a few seconds...\n"))
list_csq <- foreach(i = 1:iterations, .packages = packages) %dopar% {
  my_body <- paste0("{ \"hgvs_notations\" : [", paste0('"', hgvs_notations[(300 * (i - 1) + 1):(300 * i)], '"', collapse = ", "), " ] }")
  r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = my_body)
  fromJSON(toJSON(content(r)))
}
df_csq <- list_csq %>% rbindlist(fill = TRUE)
stopCluster(cl)

message(paste0("\nfetching completed!!!\n"))

message(paste0("\nannotation started...\n"))

# transcript annotation
annotation_transcript <- df_csq %>%
  select(input, most_severe_consequence, transcript_consequences, colocated_variants) %>%
  unnest(transcript_consequences) %>%
  select(c("input", "most_severe_consequence", "gene_id", "gene_symbol", "consequence_terms", "impact", "biotype", "transcript_id", "polyphen_prediction", "sift_prediction")) %>%
  unnest(cols = c("input", "most_severe_consequence", "gene_id", "gene_symbol", "consequence_terms", "impact", "biotype", "transcript_id", "polyphen_prediction", "sift_prediction"), keep_empty = TRUE) %>%
  filter(most_severe_consequence == consequence_terms) %>%
  select(-"consequence_terms") %>%
  group_by(input, most_severe_consequence) %>%
  summarize(transcript_id = paste(transcript_id, collapse = ","), gene_id = paste(gene_id, collapse = ","), impact = paste(impact, collapse = ","), gene_symbol = paste(gene_symbol, collapse = ","), biotype = paste(biotype, collapse = ","), polyphen_prediction = paste(polyphen_prediction, collapse = ","), sift_prediction = paste(sift_prediction, collapse = ","), .groups = "keep")


# annotation - allele frequency, clinvar significance and pubmed id
annotation_AF <- df_csq %>%
  select(input, most_severe_consequence, colocated_variants) %>%
  unnest(colocated_variants) %>%
  select(c("input", "most_severe_consequence", "id", "pubmed", "minor_allele_freq", "minor_allele", "clin_sig")) %>%
  unnest(cols = c("input", "most_severe_consequence", "id", "minor_allele_freq", "minor_allele", "clin_sig"), keep_empty = TRUE) %>%
  filter(!is.na(minor_allele_freq)) %>%
  unnest(pubmed, keep_empty = TRUE) %>%
  group_by(input, most_severe_consequence, id, minor_allele_freq, minor_allele) %>%
  summarize(clin_sig = paste(clin_sig, collapse = ","), pubmed = paste(pubmed, collapse = ","), .groups = "keep")

# merge annotation files
df_anno <- df_anno %>%
  mutate(input = paste0(`#CHROM`, ":g.", POS, REF, ">", ALT)) %>%
  left_join(annotation_AF, by = "input") %>%
  left_join(annotation_transcript, by = c("input", "most_severe_consequence")) %>%
  select(-"input")

# calculate variant allele frequency from minor allele freq
df_anno <- df_anno %>%
  mutate(variant_allele_freq.1kg = case_when(
    ALT == minor_allele ~ minor_allele_freq, 
    REF == minor_allele ~ 1 - minor_allele_freq))

# clean the final ouput
df_anno <- df_anno %>% select(c("#CHROM", "POS", "REF", "ALT", "DP", "TYPE", "RO", "AO", "AD_alt_vs_ref", "most_severe_consequence", "variant_allele_freq.1kg", "id", "minor_allele_freq", "minor_allele", "clin_sig", "pubmed","transcript_id", "gene_id", "impact", "gene_symbol", "biotype", "polyphen_prediction", "sift_prediction")) %>% rename("CHR" = "#CHROM", "total_read_depth" = "DP", "ref_read_depth" = "RO", "alt_read_depth" = "AO", "rsid" = "id", "clinvar_significance" = "clin_sig")
# write out the annotated tsv file 
df_anno %>% fwrite(output_tsv, sep = "\t")
message(paste0("\n", output_tsv," generated\n"))

# ===============================================================================
#    Extra reports
# ===============================================================================

# get summary stats   
n_variants.before_split <- my_vcf %>% nrow()
n_variants.after_split <- df_anno %>% nrow()
individual_level_columns <- c((which(colnames(my_vcf) == "FORMAT") + 1):dim(my_vcf)[2])
individual_names <- names(my_vcf)[individual_level_columns]
individual_numbers <- length(individual_names)
message(paste0("\nThere are ", n_variants.before_split, " variant sites, and ", n_variants.after_split, " variants after splitting multi-allelic variant calls into separate lines"))
message(paste0("\nThere are ", individual_numbers, " samples in the vcf: ", paste0(individual_names, collapse = ", "), "\n"))

# get metadata
if (opt$getmeta == TRUE) {
  strings <- readLines(input_vcf)
  patterns <- c("^##INFO", "^##FORMAT")
  list_metadata <- patterns %>% map(get_metadata, string = strings)
  meta_info_tsv <- paste(opt$input_dir, paste0(basename, ".info_meta.tsv"), sep = "/")
  meta_format_tsv <- paste(opt$input_dir, paste0(basename, ".format_meta.tsv"), sep = "/")
  fwrite(list_metadata[[1]], meta_info_tsv, sep = "\t")
  fwrite(list_metadata[[2]], meta_format_tsv, sep = "\t")
  message(paste0("\n", meta_info_tsv,"&", meta_format_tsv, " generated\n"))
  
}

# get all the info fields  
if (opt$getinfo == TRUE) {
  colnames(my_vcf) <- get_header(strings)
  message(paste0("\ngenerate a file containing all the variant level information... It might take a few seconds...\n"))
  fwrite(get_info(my_vcf), paste(opt$input_dir, paste0(basename, ".info.tsv"), sep = "/"), sep = "\t")
  message(paste0("\n", paste(opt$input_dir, paste0(basename, ".info.tsv"), sep = "/")," generated\n"))
  
}

# get all the FORMAT fields with different samples in separate rows
if (opt$getformat == TRUE) {
  message(paste0("\ngenerate a file containing all the individual level information...\n"))
  fwrite(get_format(my_vcf), paste(opt$input_dir, paste0(basename, ".format.tsv"), sep = "/"), sep = "\t")
  message(paste0("\n", paste(opt$input_dir, paste0(basename, ".format.tsv"), sep = "/")," generated\n"))
}

# print version information about R, the OS and attached or loaded packages for reproducibility
sessionInfo()
