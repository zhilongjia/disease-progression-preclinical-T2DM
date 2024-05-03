library(coloc)
library(dplyr)
library(readxl)
library(fs)

res_path <- 'coloc_results'
if (!dir_exists(res_path)) {
  dir_create(res_path)
}

# exposure
x_gwas_id <- read.csv('data/t2dm_gwas_files.csv')

# outcome gwas id (from opengwas)
y_gwas_id <- read.csv('data/disease_gwas_ids.csv')

pheno_x <- x_gwas_id$pheno
pheno_y <- y_gwas_id$pheno_abv
x_gwas_files <- x_gwas_id$gwas_id
pheno_y_name <- y_gwas_id$pheno

sample_size_x <- x_gwas_id$Sample_size
sample_size_y <- y_gwas_id$Sample_size
case_x <- x_gwas_id$N_case
case_y <- y_gwas_id$N_case

x_gwas_path <- 'data/gwas_summary/t2dm'
y_gwas_path <- 'data/gwas_summary/disease'

# read gwas files
i <- 1

coloc_res_path <- 'coloc_results/t2dm_disease.csv'
if (file_exists(coloc_res_path)) {
  coloc_results <- read.csv(coloc_res_path)
} else {
  coloc_results <- read.csv('coloc_results/coloc_example.csv')
}


for (i in 1:length(pheno_x)) {
  x <- pheno_x[i]
  x_gwas <- read.csv(file.path(x_gwas_path, x_gwas_files[i]), sep = ' ')
  # keep only the columns we need to save memory
  x_gwas <- x_gwas[, c('SNP', 'P', 'BETA', 'SE')]
  
  s_x <- case_x[i] / sample_size_x[i]
  
  for (j in 1:length(pheno_y)) {
    
    # check if results already exist
    coloc_xy <- coloc_results %>% filter(x == pheno_x[i] & y == pheno_y_name[j])
    if (nrow(coloc_xy) > 0) {
      print(paste0('Coloc results for ', x, ' and ', pheno_y_name[j], ' already exist. Skipping...'))
      next
    }
    
    print(paste0('Running coloc for ', x, ' and ', pheno_y_name[j]))
    y <- pheno_y[j]
    y_gwas <- read.csv(file.path(y_gwas_path, paste0(pheno_y[j], '.txt.gz')), sep = ' ')
    y_gwas <- y_gwas[, c('SNP', 'P', 'BETA', 'SE')]
    
    s_y <- case_y[j] / sample_size_y[j]
    
    input <- merge(x_gwas, y_gwas, by = 'SNP', all = FALSE, 
                   suffixes=c("_x", "_y"))
    # drop duplicate SNPs
    input <- input[!duplicated(input$SNP), ]
    result <- coloc.abf(dataset1=list(pvalues = input$P_x,
                                      type = "cc", 
                                      N = sample_size_x[i],
                                      beta = input$BETA_x,
                                      varbeta = input$SE_x^2,
                                      S = s_x, 
                                      snp = input$SNP),
                        
                        dataset2=list(pvalues=input$P_y,
                                      type="cc", 
                                      beta = input$BETA_y,
                                      varbeta = input$SE_y^2,
                                      N = sample_size_y[j],
                                      S = s_y, 
                                      snp = input$SNP))
    
    coloc_result <- result$results %>% filter(SNP.PP.H4 >= 0.75)
    # write.csv(coloc_result, 'coloc_results/coloc_example.csv', row.names = F)
    if (nrow(coloc_result) > 0) {
      print(paste0(nrow(coloc_result), ' SNPs with PP.H4 >= 0.75 for ', x, ' and ', y))
      coloc_result$x <- x
      coloc_result$y <- pheno_y_name[j]
      coloc_results <- rbind.data.frame(coloc_results, coloc_result)
      
    } else {
      # set the first row of coloc_result to NA
      coloc_result[1, ] <- NA
      coloc_result$x <- x
      coloc_result$y <- pheno_y_name[j]
      coloc_results <- rbind.data.frame(coloc_results, coloc_result)
      print(paste0('No SNPs with PP.H4 >= 0.75 for ', x, ' and ', y))
    }
    # save results
    write.csv(coloc_results, 'coloc_results/t2dm_disease.csv', row.names = F)
  }
}



