## Association Test
# Load required packages 
library(optparse)
library(data.table)
library(dplyr)
# Set seed
set.seed(10)
option_list <- list(
  make_option(c("-a", "--geno_a"), type="character", default="geno_a.geno",
    help="Insert path to vcf.gz file"),
  make_option(c("-b", "--geno_b"), type="character", default="geno_b.geno",
    help="Insert path to tsv file containing gene boundaries"),
  make_option(c("-f", "--assoc_test"), type="character", default=TRUE,
    help="Perform association tests: arguments should be TRUE or FALSE"), 
  make_option(c("-c", "--chi_file"), type="character", default="chi_result.csv",
    help="Name of file containing association test results"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",
    help="File containing metadata with DST information"),
  make_option(c("-d", "--drug"), type="character", default="rifampicin",
    help="Column name for drug DST information in metadata file")  
    )
# Parse Options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
# Perform pairwise SNP Association test
pair_snp_test <- function(geno_a, geno_b, chi_file){
  cat(" *Performing Chi-squared analysis \n")
  geno_a <- read.table(geno_a)
  geno_b <- read.table(geno_b)
  pair_tables <- vector('list', ncol(geno_a))
  for (i in 1:ncol(geno_a)){
    for (j in 1:ncol(geno_b)){
      pair_tables[[i]][[j]]<-table(geno_a[,i], geno_b[,j])
      pair_tables[[i]][[j]] <- fisher.test(pair_tables[[i]][[j]])
      pair_tables[[i]][[j]]  <- data.frame(colnames(geno_a)[i], colnames(geno_b)[j], pair_tables[[i]][[j]][[1]], pair_tables[[i]][[j]][[2]][1],pair_tables[[i]][[j]][[2]][2],pair_tables[[i]][[j]][[3]])
      colnames(pair_tables[[i]][[j]])<- c("snp1","snp2","pair_pvalue","pair_lowerCI", "pair_upperCI","pair_OR")
    }
  }
  for (i in 1:length(pair_tables)){
    pair_tables[[i]] <- do.call("rbind", pair_tables[[i]])
  }
  pair_tables <- do.call("rbind", pair_tables)
  pair_tables$p.adj <- p.adjust(as.numeric(as.character(pair_tables$pair_pvalue)), method="fdr")
  return(pair_tables)
}
# Function to perform SNP:phenotype Association test
snp_pheno_test <- function(geno, meta, drug){
  metadata <- read.csv(meta, header=TRUE)
  metadata <- metadata[!is.na(metadata[,drug]),]
  geno <- read.table(geno)
  geno <- geno[rownames(geno) %in% metadata$id,]
  phen_tables <- list()
  for (i in 1:ncol(geno)){
    phen_tables[[i]]<-table(geno[,i], metadata[,drug])
    if (nrow(phen_tables[[i]]) >1 & ncol(phen_tables[[i]])>1){
      phen_tables[[i]]<- fisher.test(phen_tables[[i]])
      phen_tables[[i]] <- data.frame(colnames(geno)[i], phen_tables[[i]][[1]], phen_tables[[i]][[2]][1],phen_tables[[i]][[2]][2],phen_tables[[i]][[3]])
      colnames(phen_tables[[i]])<- c("snp","pvalue","lowerCI", "upperCI","OR")
    }
    else{
    phen_tables[[i]]<- data.frame("NA","NA","NA","NA","NA")
    colnames(phen_tables[[i]])<- c("snp","pvalue","lowerCI", "upperCI","OR")
    }
  }
  phen_tables <- rbindlist(phen_tables)
  phen_tables$p.adj <- p.adjust(as.numeric(as.character(phen_tables$pvalue)), method="fdr")
  return(phen_tables)
}
# Function to combine Chi-squared results into one file
combine_results <- function(tables, phen_tables_a, phen_tables_b, chi_filename){
  colnames(phen_tables_a) <- c("snp1","snp1_pvalue","snp1_lowerCI", "snp1_upperCI","snp1_OR", "snp1_padj")
  colnames(phen_tables_b) <- c("snp2","snp2_pvalue","snp2_lowerCI", "snp2_upperCI","snp2_OR", "snp2_padj")
  chi <- full_join(tables, phen_tables_a, by="snp1")
  chi <- full_join(chi, phen_tables_b, by="snp2")
  write.csv(chi, chi_filename, row.names=FALSE)
}
# Run
if(opt$assoc_test==TRUE){
  tables <- pair_snp_test(opt$geno_a, opt$geno_b, chi_file)
  phen_tables_a <- snp_pheno_test(opt$geno_a, opt$meta, opt$drug)
  phen_tables_b <- snp_pheno_test(opt$geno_b, opt$meta, opt$drug)
  combine_results(tables, phen_tables_a, phen_tables_b, opt$chi_file)
}
