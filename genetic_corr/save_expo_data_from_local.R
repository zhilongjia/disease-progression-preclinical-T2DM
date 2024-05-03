rm(list = ls())
library(TwoSampleMR)
library(ieugwasr)
library(fs)

source('ld_clump_local.R')

expo_gwas_id <- read.csv('data/FAG.csv')

expos <- expo_gwas_id$pheno
filenames <- expo_gwas_id$pheno_abv

save_path <- 'exposure/FAG'
gwas_path <- 'data/gwas_summary/FAG/tophits/'
# i <- 10

if (!dir_exists(save_path)) {
  dir_create(save_path)
}

for (i in 1:length(expos)) {
  print(sprintf('extract for %s', expos[i]))
  filename <- sprintf('%s/%s_top.csv', gwas_path, filenames[i])
  # read exposure data

  if (expos[i] == 'T2DM' || expos[i] == 'T1DM') {
    filename <- sprintf('%s/%s_top.txt', gwas_path, filenames[i])
    expo <- read_exposure_data(
      filename = filename,
      sep = ' ',
      clump = F,
      snp_col = 'SNP',
      chr_col = 'CHR',
      se_col = 'SE',
      beta_col = 'BETA',
      effect_allele_col = 'A1',
      other_allele_col = 'A2',
      pval_col = 'P',
      eaf_col = 'EAF'
    )
  } else {
    expo <- read_exposure_data(
      filename = filename,
      sep = ',',
      clump = F,
      snp_col = 'SNP',
      chr_col = 'CHR',
      se_col = 'SE',
      beta_col = 'BETA',
      effect_allele_col = 'A1',
      other_allele_col = 'A2',
      pval_col = 'P',
      eaf_col = 'EAF'
    )
  }
  
  # expo$id.exposure <- 'Anorexia_nervosa'
  expo$id.exposure <- expos[i]
  
  # f statistic > 10
  expo$fstat <- (expo$beta.exposure / expo$se.exposure)^2
  expo <- subset.data.frame(expo, fstat >= 10)
  
  # LD clump using EUR reference panel
  expo$rsid <- expo$SNP
  expo$pval <- expo$pval.exposure
  expo_clumped <- ld_clump_local(expo, clump_kb = 1000,
                                 clump_r2 = 0.001, clump_p = 1,
                                 plink_bin = 'reference/plink',
                                 bfile = 'reference/EUR')
  # expo_clumped <- clump_data(expo, clump_kb = 1000, clump_r2 = 0.1)
  
  expo_filename <- sprintf('%s/%s_expo.csv', save_path, expos[i])
  write.csv(expo_clumped, expo_filename, row.names = F)
  
}



