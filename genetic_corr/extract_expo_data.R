library(TwoSampleMR)
library(ieugwasr)
library(fs)

source('ld_clump_local.R')

save_path <- 'exposure/t2dm'

# make sure the directory exists
if (!dir_exists(save_path)) {
  dir_create(save_path)
}

# read exposure data
filename <- sprintf('data/gwas_summary/t2dm/res_subtype2_pc10_top.csv')

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
  eaf_col = 'MAF',
)

expo$id.exposure <- 'Subtype2'

# f statistic > 10
expo$fstat <- (expo$beta.exposure / expo$se.exposure)^2
expo <- subset.data.frame(expo, fstat >= 10)

# LD clump using EUR reference panel
expo$rsid <- expo$SNP
expo$pval <- expo$pval.exposure
expo_clumped <- ld_clump_local(expo, clump_kb = 1000,
                               clump_r2 = 0.1, clump_p = 1,
                               plink_bin = 'reference/plink',
                               bfile = 'reference/EUR')

# expo_clumped <- clump_data(expo, clump_kb = 1000, clump_r2 = 0.1)
print(length(expo_clumped$SNP))

expo_filename <- sprintf('%s/%s_expo.csv', save_path, 'T2DM_Subtype2')
write.csv(expo_clumped, expo_filename, row.names = F)


