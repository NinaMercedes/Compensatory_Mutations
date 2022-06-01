## Pre-process Data
# Load required packages 
library(optparse)
library(data.table)
library(dplyr)
# Set seed
set.seed(10)
option_list <- list(
  make_option(c("-v", "--vcf_gz"), type="character", default="geno_a.vcf.gz",
    help="Insert path to vcf.gz file"),
  make_option(c("-r", "--region"), type="character", default="region.tsv",
    help="Insert path to tsv file containing gene boundaries"),
  make_option(c("-f", "--filter_vcf"), type="character", default=TRUE,
    help="Filter vcf file: arguments should be TRUE or FALSE"),
  make_option(c("-o", "--out_vcf"), type="character", default=FALSE,
    help="Name of filtered vcf.gz to output for gene "), 
  make_option(c("-j", "--make_geno"), type="character", default=TRUE,
    help="Make genotype file: arguments should be TRUE or FALSE"), 
  make_option(c("-g", "--out_geno"), type="character", default="geno_a.geno",
    help="Name of geno file to output"),
  make_option(c("-s", "--sample_names"), type="character", default="sample_names.txt",
    help="Name of file containing sample names"), 
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",
    help="File containing metadata with DST information"),
  make_option(c("-d", "--drug"), type="character", default="rifampicin",
    help="Column name for drug DST information in metadata file"), 
  make_option(c("-l", "--rem_lineage"), type="character", default=TRUE,
    help="Remove lineage SNPs: arguments should be TRUE or FALSE"), 
  make_option(c("-n", "--out_geno_nolin"), type="character", default="geno_a_nolin.geno",
    help="Name of geno file to output with no lineage SNPs"),
  make_option(c("-x", "--missense"), type="character", default=TRUE,
    help="Filter missense SNPs: arguments should be TRUE or FALSE"))
# Parse Options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
# Function to filter vcf by gene region and get missense SNPs
filter_vcf <- function(region, vcf_gz, out_vcf, missense){
  cat(" *Filtering vcf \n")
  if(missense==TRUE){
  system(paste("bcftools view -Oz -R", region, "-i 'INFO/BCSQ[*] ~ \"missense*\"'", vcf_gz, " > temp.vcf.gz"))
  }
  else{
  system(paste("bcftools view -Oz -R", region, vcf_gz, " > temp.vcf.gz"))
  }
  system(paste("bcftools index -c --force temp.vcf.gz"))
  system(paste("bcftools norm -Oz -m-any temp.vcf.gz >", out_vcf))
  system(paste("rm temp.vcf.gz"))
}
# Function to make genotype file
make_geno <- function(out_vcf, geno, sample_names){
  cat(" *Making geno file  \n")
  system(paste("bcftools index -c --force",out_vcf))
  system(paste("bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]\n'", out_vcf, " | tr '|' '/' | sed 's/\\.\\/\\./Ns/g' | sed 's/0\\/1/0.5/g' | sed 's/[123456789]\\/[123456789]/1/g' | sed 's/0\\/0/0/g' >", geno))
  system(paste("bcftools query -l", out_vcf, " >", sample_names))
  system(paste("awk \'BEGIN{OFS=\"_\"} {print $1,$2,$3}\'", geno, "> tmp"))
  system(paste("awk \'NR==FNR{a[NR]=$0;next}{$1=a[FNR]}1\' tmp ", geno, ">tmp2 &&  cut --complement -d' ' -f 2,3 tmp2 >tmp3 && mv tmp3", geno, "&& rm tmp2 && rm tmp"))  
}
# Function to format genotype
format_geno <- function(geno_table, sample_names){
  cat(" *Formatting genotypic data \n")
  geno_table <- read.table(geno_table)
  sample_name <- read.table(sample_names)
  snp_name <- geno_table[,1]
  geno_table <- geno_table[,-c(1)]
  is.na(geno_table)<- geno_table == "N" 
  geno_table <- t(geno_table)
  colnames(geno_table)<- snp_name
  geno_table <- apply(geno_table,2, function(x) as.numeric(as.character(x)))
  geno_table <- geno_table[,-(which(colSums(geno_table)==0))]
  geno_table <- data.frame(geno_table)
  rownames(geno_table)<- sample_name$V1
  return(geno_table)
}
# Functions to remove lineage SNPs
lineage_loop <- function(geno_table, metadata, drug){
  metadata <- metadata[grep("1", metadata[,drug]),]
  lineage <- matrix(ncol=1, nrow=ncol(geno_table))
  lineage <- data.frame(lineage)
  mylist <- list()
  for (i in colnames(geno_table)){
    mylist[[i]] <- geno_table[grep("1",geno_table[,i]),]
    mylist[[i]] <- rownames(mylist[[i]])
    mylist[[i]] <- metadata[metadata$id %in% mylist[[i]],]
    mylist[[i]]<- nrow(unique(mylist[[i]]['lineage']))
    if (mylist[i] < 2 ){lineage[i,1] <- names(mylist[i])
    } 
  }
return(lineage)
}
rem_lin <- function(metadata, drug, geno, geno_nolin_file, sample_names){
  metadata <- read.csv(metadata, header=TRUE)
  metadata$lineage <- gsub("lineage1.*", "lineage1", metadata$lineage)
  metadata$lineage <- gsub("lineage2.*", "lineage2", metadata$lineage)
  metadata$lineage <- gsub("lineage3.*", "lineage3", metadata$lineage)
  metadata$lineage <- gsub("lineage4.*", "lineage4", metadata$lineage)
  metadata$lineage <- gsub("lineage5.*", "lineage5", metadata$lineage)
  metadata$lineage <- gsub("lineage6.*", "lineage6", metadata$lineage)
  metadata$lineage <- gsub("lineage7.*", "lineage7", metadata$lineage)
  geno <- format_geno(geno, sample_names)
  lineage <- lineage_loop(geno, metadata, drug)
  lineage <- data.frame(lineage)
  colnames(lineage) <- c("lin")
  lineage <- lineage[!is.na(lineage$lin),]
  geno_no_lin <- geno[,-which(colnames(geno) %in% lineage)]
  write.table(geno_no_lin, geno_nolin_file, sep = " ", quote=FALSE, row.names=TRUE)
}
# Run
if(opt$filter_vcf==TRUE){
  filter_vcf(opt$region, opt$vcf_gz, opt$out_vcf, opt$missense)
}
if(opt$make_geno==TRUE){
  make_geno(opt$out_vcf, opt$out_geno, opt$sample_names)
}
if(opt$rem_lin==TRUE){
  rem_lin(opt$metadata, opt$drug, opt$out_geno, opt$out_geno_nolin, opt$sample_names)
}