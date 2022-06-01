# Interpret Results
# Load required packages 
library(optparse)
library(data.table)
library(dplyr)
# Set seed
set.seed(10)
option_list <- list(
  make_option(c("-i", "--interpret_results"), type="character", default=TRUE,
    help="Interpret results: Arguments TRUE or FALSE"),
  make_option(c("-r", "--med_file"), type="character", default="med_result.csv",
    help="Name of file containing mediation and association test results"),
  make_option(c("-a", "--info_file_r1"), type="character", default="info1.csv",
    help="Name of file containing region 1 information"),
  make_option(c("-b", "--info_file_r2"), type="character", default="info2.csv",
    help="Name of file containing region 2 information"),
  make_option(c("-o", "--out_file"), type="character", default="info.csv",
    help="Name of file to output significant results"),
  make_option(c("-c", "--chi_only"), type="character", default="chi_only.csv",
    help="Interpret results: Arguments TRUE or FALSE")
    )
# Parse Options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
# Function to interpret and filter results
interpret_results <- function(results_file, info_file_r1, info_file_r2, out_file, chi_file){
  #read in data
  results <- read.csv(results_file, header=TRUE) 
  info_r1 <- read.csv(info_file_r1, header=TRUE)
  info_r2 <- read.csv(info_file_r2, header=TRUE)
  #join together results and snp information- see flow diagram
  results[,-c(1:3)] <- apply(results[,-c(1:3)], 2, function(x) as.numeric(as.character(x)))
  results <- left_join(results, info_r1)
  results <- left_join(results, info_r2)
  #Get significant mediation results and categorise
  results_sig <- results %>% filter((ab_ci_l>0 & ab_ci_u>0)|(ab_ci_l<0 & ab_ci_u<0))
  results_sig$mediation_type <- ifelse(results_sig$c_ci_l<0 & results_sig$c_ci_u>0, "Full", "Partial")
  results_sig$mediation_type2 <- ifelse((results_sig$mediation_type=="Partial" & results_sig$ab_est>0 & results_sig$c_est>0)|(results_sig$mediation_type=="Partial" & results_sig$ab_est<0 & results_sig$c_est<0),"Complementary", "Competitive")
  results_sig$ab_tot <- results_sig$ab_std/results_sig$tot_std
  results_sig$significant_analyses <- ifelse(results_sig$p.adj <0.05 & results_sig$snp1_padj <0.05  & results_sig$snp2_padj <0.05, "Association and Mediation", "Mediation")
  results_sig$putative_compensatory <- ifelse(results_sig$mediation_type2=="Complementary" & results_sig$ab_tot>0.05, "Highly likely", "Likely")
  results_sig <- within(results_sig, mediation_type2[mediation_type=="Full"] <- "NA")
  results_sig <- within(results_sig, putative_compensatory[mediation_type2=="Competitive"] <- "Unlikely")
  results_sig <- within(results_sig, putative_compensatory[mediation_type2=="Complementary" & significant_analyses=="Mediation"] <- "Highly unlikely")
  results_sig <- within(results_sig, putative_compensatory[mediation_type2=="Competitive" & significant_analyses=="Mediation"] <- "Highly unlikely")
  results_sig <- within(results_sig, putative_compensatory[mediation_type=="Full" & ab_tot>0.05] <- "Highly likely")
  results_sig <- within(results_sig, putative_compensatory[mediation_type=="Full" & ab_tot<0.05] <- "Unlikely")
  write.csv(results_sig, out_file, row.names=FALSE)
  #Get significant results for association test only
  results_notsig <- results %>% filter((ab_ci_l<0 & ab_ci_u>0))
  chi_only <- results_notsig %>% filter(results_notsig$p.adj <0.05 & results_notsig$snp1_padj <0.05  & results_notsig$snp2_padj <0.05)
  chi_only$significant_analyses <- "Association"
  write.csv(chi_only, chi_file, row.names=FALSE)
}
# Run
if(opt$interpret_results==TRUE){
  interpret_results(opt$med_file, opt$info_file_r1, opt$info_file_r2, opt$out_file, opt$chi_only)
}
