## Mediation Test
#Load required packages 
library(optparse)
library(data.table)
library(dplyr)
library(lavaan)
# Set seed
set.seed(10)
option_list <- list(
  make_option(c("-a", "--geno_a"), type="character", default="geno_a.geno",
    help="Insert path to vcf.gz file"),
  make_option(c("-b", "--geno_b"), type="character", default="geno_b.geno",
    help="Insert path to tsv file containing gene boundaries"),
  make_option(c("-t", "--med_test"), type="character", default=TRUE,
    help="Perform association tests: arguments should be TRUE or FALSE"), 
  make_option(c("-s", "--sample_names"), type="character", default="sample_names.txt",
    help="Name of file containing sample names"),
  make_option(c("-c", "--chi_file"), type="character", default="chi_result.csv",
    help="Name of file containing association test results"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",
    help="File containing metadata with DST information"),
  make_option(c("-d", "--drug"), type="character", default="rifampicin",
    help="Column name for drug DST information in metadata file"),
  make_option(c("-r", "--med_file"), type="character", default="med_result.csv",
    help="Name of file containing mediation test results")
    )
# Parse Options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
# Function to perform mediation analysis
mediation <- function(geno_a, geno_b, chifile, metadata, mediation_file, sample_names, phenotype){
  geno_a <- read.table(geno_a)
  geno_b <- read.table(geno_b)
  chi <- read.csv(chifile, header=TRUE)
  metadata <- read.csv(metadata, header=TRUE)
  model_summary <- list()
  chi$a_est <- NA
  chi$a_se <- NA
  chi$a_z <- NA
  chi$a_p <- NA
  chi$a_ci_l <- NA
  chi$a_ci_u <- NA
  chi$a_std <- NA
  chi$b_est <- NA
  chi$b_se <- NA
  chi$b_z <- NA
  chi$b_p <- NA
  chi$b_ci_l <- NA
  chi$b_ci_u <- NA
  chi$b_std <- NA
  chi$c_est <- NA
  chi$c_se <- NA
  chi$c_z <- NA
  chi$c_p <- NA
  chi$c_ci_l <- NA
  chi$c_ci_u <- NA
  chi$c_std <- NA
  chi$ab_est <- NA
  chi$ab_se <- NA
  chi$ab_z <- NA
  chi$ab_p <- NA
  chi$ab_ci_l <- NA
  chi$ab_ci_u <- NA
  chi$ab_std <- NA
  chi$tot_est <- NA
  chi$tot_se <- NA
  chi$tot_z <- NA
  chi$tot_p <- NA
  chi$tot_ci_l <- NA
  chi$tot_ci_u <- NA
  chi$tot_std <- NA
  for (i in 1:nrow(chi)){
    snp1 = as.character(chi[i,"snp1"])
    print(snp1)
    snp2 = as.character(chi[i,"snp2"])
    print(snp2)
    tg <- data.frame(geno_a[,snp1] ,geno_b[,snp2])
    sample_name <- read.table(sample_names, header=FALSE)
    tg$id <-sample_name$V1
    drug <- data.frame(metadata$id, metadata[,phenotype])
    colnames(drug) <- c("id","drug")
    tg <- left_join(tg, drug)
    colnames(tg) <- c("snp1", "snp2", "id", "drug")
    tg <- tg[!is.na(tg$drug),]
    if(length(unique(tg$drug))>=2 & length(unique(tg$snp1))>=2 & length(unique(tg$snp2))>=2){
      if (all(tg$snp1 == tg$snp2)==FALSE){
        model = '
        #direct effect
        drug ~ c*snp1
        #mediator
        snp2 ~ a*snp1
        drug ~ b*snp2
        #indirect effect
        ab := a*b
        #total effect
        total := c + (a*b)
        '
        nBoots <- 1000
        fit <- sem(model, data = tg, missing = 'listwise', se = "boot", estimator='dwls', link='probit', bootstrap = nBoots)
        if (fit@Fit@converged==TRUE){
          model_summary[[i]] = summary(fit,fit.measures=TRUE, rsquare=TRUE, standardized=TRUE, ci=TRUE)
          chi[i,]$c_est <- model_summary[[i]]$PE[1,6]
          chi[i,]$c_se <- model_summary[[i]]$PE[1,7]
          chi[i,]$c_z <- model_summary[[i]]$PE[1,8]
          chi[i,]$c_p <- model_summary[[i]]$PE[1,9]
          chi[i,]$c_ci_l <- model_summary[[i]]$PE[1,10]
          chi[i,]$c_ci_u <- model_summary[[i]]$PE[1,11]
          chi[i,]$c_std <- model_summary[[i]]$PE[1,13]
          chi[i,]$a_est <- model_summary[[i]]$PE[2,6]
          chi[i,]$a_se <- model_summary[[i]]$PE[2,7]
          chi[i,]$a_z <- model_summary[[i]]$PE[2,8]
          chi[i,]$a_p <- model_summary[[i]]$PE[2,9]
          chi[i,]$a_ci_l <- model_summary[[i]]$PE[2,10]
          chi[i,]$a_ci_u <- model_summary[[i]]$PE[2,11]
          chi[i,]$a_std <- model_summary[[i]]$PE[2,13]
          chi[i,]$b_est <- model_summary[[i]]$PE[3,6]
          chi[i,]$b_se <- model_summary[[i]]$PE[3,7]
          chi[i,]$b_z <- model_summary[[i]]$PE[3,8]
          chi[i,]$b_p <- model_summary[[i]]$PE[3,9]
          chi[i,]$b_ci_l <- model_summary[[i]]$PE[3,10]
          chi[i,]$b_ci_u <- model_summary[[i]]$PE[3,11]
          chi[i,]$b_std <- model_summary[[i]]$PE[3,13]
          chi[i,]$ab_est <- model_summary[[i]]$PE[7,6]
          chi[i,]$ab_se <- model_summary[[i]]$PE[7,7]
          chi[i,]$ab_z <- model_summary[[i]]$PE[7,8]
          chi[i,]$ab_p <- model_summary[[i]]$PE[7,9]
          chi[i,]$ab_ci_l <- model_summary[[i]]$PE[7,10]
          chi[i,]$ab_ci_u <- model_summary[[i]]$PE[7,11]
          chi[i,]$ab_std <- model_summary[[i]]$PE[7,13]
          chi[i,]$tot_est <- model_summary[[i]]$PE[8,6]
          chi[i,]$tot_se <- model_summary[[i]]$PE[8,7]
          chi[i,]$tot_z <- model_summary[[i]]$PE[8,8]
          chi[i,]$tot_p <- model_summary[[i]]$PE[8,9]
          chi[i,]$tot_ci_l <- model_summary[[i]]$PE[8,10]
          chi[i,]$tot_ci_u <- model_summary[[i]]$PE[8,11]
          chi[i,]$tot_std <- model_summary[[i]]$PE[8,13]
        }
      }
    }
  }
  write.csv(chi, mediation_file)
}
#Run
if(opt$med_test==TRUE){
  mediation(opt$geno_a, opt$geno_b, opt$chi_file, opt$metadata, opt$med_file, opt$sample_names, opt$drug)
}