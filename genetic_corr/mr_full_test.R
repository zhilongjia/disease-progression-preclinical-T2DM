library(TwoSampleMR)
library(ieugwasr)
library(MRPRESSO)

source('ld_clump_local.R')
source('mr_modified.R')

# load mr_example.RData
# load('mr_example.RData')

mr_full_test <- function(dat, expo, outcome, presso_test = T) {
  if ((nrow(dat) == 0) || nrow(dat[dat$mr_keep == T, ] ) == 0) {
    # no snp present for MR analysis
    print('Not enough snp can be analysed.')
    return(NULL)
  } else {
    dat <- dat[dat$mr_keep == T, ]
    if (nrow(dat) == 1) {
      mr_method <- c('mr_wald_ratio')
    } else if (nrow(dat) == 2) {
      mr_method <- c('mr_ivw')
    } else {
      mr_method <- c('mr_ivw', 'mr_egger_regression', 'mr_raps', 
                     'mr_two_sample_ml',
                     'mr_weighted_median', 'mr_weighted_mode')
    }
  }
  
  mr_res <- mr_modified(dat, method_list = mr_method)
  
  if (nrow(dat) >= 3) {
    res_pleiotropy <- mr_pleiotropy_test(dat)
    res_heterogeneity <- mr_heterogeneity(dat, method_list = 'mr_ivw')
    if (nrow(dat) >= 4 && presso_test == T) {
      # try mr presso
      res_presso <- try(mr_presso('beta.outcome', 'beta.exposure', 'se.outcome',
                                  'se.exposure', data = dat))
      # if mr presso failed, set to NULL
      if ("try-error" %in% class(res_presso)) {
        res_presso_mr <- NULL
        res_presso_pleiotropy <- NULL
      } else {
        # save pleiotropy test results and presso results
        res_presso_mr <- res_presso$`Main MR results`
        res_presso_mr$id.exposure <- expo
        res_presso_mr$id.outcome <- outcome

        res_presso_pleiotropy <- data.frame(res_presso$`MR-PRESSO results`$`Global Test`)
        res_presso_pleiotropy$id.exposure <- expo
        res_presso_pleiotropy$id.outcome <- outcome
      }
      
    } else {
      res_presso_pleiotropy <- NULL
      res_presso_mr <- NULL
    }
    
  } else {
    # pleiotropy test, heterogeneity test, mr presso not performed
    res_pleiotropy <- NULL
    res_heterogeneity <- NULL
    res_presso_pleiotropy <- NULL
    res_presso_mr <- NULL
  }
  
  return(list(res = mr_res, res_pleiotropy = res_pleiotropy, 
              res_heterogeneity = res_heterogeneity,
              res_presso_mr = res_presso_mr,
              res_presso_pleiotropy = res_presso_pleiotropy))
}


# save results
# write.csv(res, 'results/mr_res_example.csv', row.names = F)
# write.csv(res_pleiotropy, 'results/pleiotropy_res_example.csv', row.names = F)
# write.csv(res_heterogeneity, 'results/heterogeneity_res_example.csv', row.names = F)
# write.csv(res_presso_mr, 'results/presso_mr_res_example.csv', row.names = F)
# write.csv(res_presso_pleiotropy, 'results/presso_pleiotropy_res_example.csv', row.names = F)
# 
# mr_method_lists <- mr_method_list()
# print(res_presso)

