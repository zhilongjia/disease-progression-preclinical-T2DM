rm(list = ls())
library(TwoSampleMR)
library(ieugwasr)
library(readr)

source('mr_full_test.R')

# exposure
expo_gwas_id <- read.csv('data/t2dm_gwas_files.csv')

# outcome gwas id (from opengwas)
out_gwas_id <- read.csv('data/disease_gwas_ids.csv')

# pheno names 
expos <- expo_gwas_id$pheno
expo_filenames <- expo_gwas_id$pheno
outcomes <- out_gwas_id$pheno
bd_names <- out_gwas_id$pheno_abv

res_path <- 'results/t2dm_disease'

# create result path if needed
if (!dir.exists(res_path)) {
  dir.create(res_path)
}

# results path for mr analysis
mr_res_path <- sprintf('%s/mr_res.csv', res_path) 
pleiotropy_res_path <- sprintf('%s/pleiotropy_res.csv', res_path)
heterogeneity_res_path <- sprintf('%s/heterogeneity_res.csv', res_path)
# presso_mr_res_path <- sprintf('%s/presso_mr_res.csv', res_path)
# presso_pleiotropy_res_path <- sprintf('%s/presso_pleiotropy_res.csv', res_path)

# load results if exists
if (file.exists(mr_res_path)) {
  res_all <- read.csv(mr_res_path)
  res_pleiotropy <- read.csv(pleiotropy_res_path)
  res_heterogeneity <- read.csv(heterogeneity_res_path)
  # res_presso_mr <- read.csv(presso_mr_res_path)
  # res_presso_pleiotropy <- read.csv(presso_pleiotropy_res_path)
  
} else {
  # load example file
  res_all <- read.csv('results/mr_example/mr_res_example.csv')
  res_pleiotropy <- read.csv('results/mr_example/pleiotropy_res_example.csv')
  res_heterogeneity <- read.csv('results/mr_example/heterogeneity_res_example.csv')
  # res_presso_mr <- read.csv('results/mr_example/presso_mr_res_example.csv')
  # res_presso_pleiotropy <- read.csv('results/mr_example/presso_pleiotropy_res_example.csv')
  
}


# diabetes
for (i in 1:length(expos)) {
  # extract exposure data, warning, allele as T should not be read as TRUE
  expo_clamped <- read.csv(sprintf('exposure/t2dm_rm_confounder/%s_expo_rm_confounder.csv',
                                   expo_filenames[i]), tryLogical = F)
  
  # outcome data: phenotypes

  j <- 22
  # diseases
  for (j in 1:length(outcomes)) {
    print(sprintf('%d: analying %s on %s', i, expos[i], outcomes[j]))
    bd <- bd_names[j]
    
    # check if the analysis has been done
    x <- subset.data.frame(res_all, id.exposure == expos[i] & id.outcome == outcomes[j])
    if (nrow(x) > 0) {
      print(sprintf('Done. Analying %s on %s', expos[i], outcomes[j]))
      next
    }
    
    # GWAS summary data for outcome
    pheno <- try(read_outcome_data(
      filename = sprintf('data/gwas_summary/disease/%s.txt.gz', bd),
      snps = expo_clamped$SNP,
      snp_col = 'SNP',
      sep = ' ',
      phenotype_col = '',
      beta_col = 'BETA',
      se_col = 'SE',
      pval_col = 'P',
      effect_allele_col = 'A1',
      other_allele_col = 'A2',
      chr_col = 'CHR',
      pos_col = 'BP',
      eaf_col = 'EAF',
      log_pval = F,
    ))
    
    
    # No. SNPs find in both exposure data and outcome data
    if (is.null(pheno)|| ("try-error" %in% class(pheno))) {
      print(sprintf('No snp matches for exposure %s and outcome %s', expos[i], outcomes[j]))
      next
    }
    
    # set outcome name
    pheno$id.outcome <- outcomes[j]
    
    # harmonise data, 
    dat <- harmonise_data(
      exposure_dat = expo_clamped,
      outcome_dat = pheno,
      action = 2
    )
    
    # mr analysis
    mr_analysis <- mr_full_test(dat, expos[i], outcomes[j], presso_test = F)
    
    # save results
    if (!is.null(mr_analysis)) {
      # load mr analysis results
      res_mr <- mr_analysis$res
      res_plei <- mr_analysis$res_pleiotropy
      res_het <- mr_analysis$res_heterogeneity
      # res_psmr <- mr_analysis$res_presso_mr
      # res_ps_plei <- mr_analysis$res_presso_pleiotropy
      
      # combine results
      res_all <- rbind(res_all, res_mr)
      if (!is.null(res_plei)) {
        res_pleiotropy <- rbind(res_pleiotropy, res_plei)
      }
      if (!is.null(res_het)) {
        res_heterogeneity <- rbind(res_heterogeneity, res_het)
      }
      # if (!is.null(res_psmr)) {
      #   res_presso_mr <- rbind(res_presso_mr, res_psmr)
      # }
      # if (!is.null(res_ps_plei)) {
      #   res_presso_pleiotropy <- rbind(res_presso_pleiotropy, res_ps_plei)
      # }
    }
    
    # save results to csv file
    write.csv(res_all, mr_res_path, row.names = F)
    write.csv(res_pleiotropy, pleiotropy_res_path, row.names = F)
    write.csv(res_heterogeneity, heterogeneity_res_path, row.names = F)
    # write.csv(res_presso_mr, presso_mr_res_path, row.names = F)
    # write.csv(res_presso_pleiotropy, presso_pleiotropy_res_path, row.names = F)
    
  }
}



