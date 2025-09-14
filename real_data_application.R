array.id <- 1
# array id: 0-7
# 0-3 for pheno LDL
# 4-7 for pheno BMI
# 0,4: snpboost with fixed step length 0.1
# 1,5: snpboost with adaptive step length
# 2,6: snpboostlss with fixed step length 0.1
# 3,7: snpboostlss with adaptive step length

library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)
library(readr)
library(gamlss)
library(pgenlibr)

setting <- data.frame(array_id = array.id)
setting <- setting %>%
  mutate(
    pheno = ifelse(array_id<4, "LDL", "BMI"),
    method = case_when(
      array_id %% 4 == 0 ~ "snpboost_fsl",
      array_id %% 4 == 1 ~ "snpboost_asl",
      array_id %% 4 == 2 ~ "snpboostlss_fsl",
      TRUE ~ "snpboostlss_asl",
    )
  )

if((setting$array_id %% 4) < 2){
  # load functions for snpboost from https://github.com/boost-PRS/snpboost
  source("../snpboost/lm_boost.R")
  source("../snpboost/snpboost.R")
  source("../snpboost/functions_snpboost.R")
  source("../snpboost/functions_snpnet.R")
}else{
  # load functions for snpboostlss from https://github.com/boost-PRS/snpboostlss
  source("../snpboostlss/lmlss_boost.R")
  source("../snpboostlss/snpboostlss.R")
  source("../snpboostlss/functions_snpboostlss.R")
  source("../snpboostlss/functions_snpboost.R")
  source("../snpboostlss/functions_snpnet.R")
}

# directory paths
dir.real.data <- "../real_data"
# absolute path of working directory for current real data analysis
dir.pheno <- paste0(dir.real.data, "/", setting$pheno)
if (!file.exists(dir.pheno)){
  dir.create(dir.pheno)
}
dir <- paste0(dir.pheno, "/", setting$method)
if (!file.exists(dir)){
  dir.create(dir)
}

# set directory to store results
dir.res <- paste0(dir.real.data, "/results")
if (!file.exists(dir.res)){
  dir.create(dir.res)
}

setwd(dir.real.data)

###  Fitting process starts ###
genotype.pfile <- "ukb_white_british_bl"  # genotype file
phenotype.file <- "ukb_white_british_bl.phe" # phenotype file, .phe, tabular data
phenotype <- setting$pheno # name of phenotype column in phenotype file

covariates <- NULL # name of covariates for mu
covariates.mu <- NULL # name of covariates for mu
covariates.sigma <- NULL # name of covariates for sigma

configs <- list(plink2.path = "../plink2/plink2", # path to plink file
                zcat.path = "cat",
                results.dir = dir, # results folder
                save = FALSE, # if true, results are saved after each batch
                prevIter = 0,
                missing.rate = 0.1,
                MAF.thresh = 0.01,
                num.snps.batch = 1000,  # not important here
                early.stopping = TRUE, # must be true to stop with validation set
                stopping.lag = 2, # not important here
                verbose = FALSE,
                mem=16000,
                niter=30000, # not important here
                nCores=16,
                standardize.variant=TRUE,
                vzs=FALSE # avoid zstd
)

start.time <- Sys.time()

if((setting$array_id %% 4) < 2) {
  fit <- snpboost(
    genotype.pfile = genotype.pfile,
    phenotype.file = phenotype.file,
    phenotype = phenotype,
    covariates = covariates,
    configs = configs,
    split.col = "split", ### name of column in phenotype file indication train, val, test set
    p_batch = 1000, # batch size
    m_batch = 1000, # maximum number of boosting steps per batch
    b_max = 20000, # maximum number of batches
    b_stop = 2, # stop if performance not improving for more than 2 batches
    sl = if(setting$array_id%%4==0){0.1}else{NULL}, # learning rate (default=0.1)
    give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
    family = "gaussian", # family argument: currently only gaussian
    metric = "MSEP" # loss function/metric: currently only negative log-likelihood for gaussian
  )
} else {
  fit <- snpboostlss(
    data.dir = paste0(dir.real.data, "/"),
    genotype.pfile = genotype.pfile,
    phenotype.file = phenotype.file,
    phenotype = phenotype,
    covariates.mu = covariates.mu,
    covariates.sigma = covariates.sigma,
    configs = configs,
    split.col = "split", ### name of column in phenotype file indication train, val, test set (named "train", "val", "test")
    p_batch = 1000, # batch size
    m_batch = 1000, # maximum number of boosting steps per batch
    b_max = 20000, # maximum number of batches
    b_stop = 2, # stop if performance not improving for more than 2 batches
    sl= if(setting$array_id%%4==2){0.1}else{NULL}, # learning rate (default=NULL, i.e. adaptive sl).
    # A value different from NULL means using the same fixed sl for mu and sigma.
    give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
    family = "gaussian", # family argument: currently only gaussian
    metric = "loss" # loss function/metric: currently only negative log-likelihood for gaussian
  )
}

end.time <- Sys.time()

run.time <- difftime(end.time, start.time, unit="mins")

paste(setting$method, "DONE. It takes", round(run.time,4), units(run.time))


### Performance measure calculation starts ###
# calculate performance measures
res <- list()
## number of included variants
coef.mu <- get_coefficients(
  fit$beta,
  fit$n_step,
  covariates=fit$configs[['covariates']]
)
res$coef.mu <- coef.mu
res$sparsity.mu <- length(coef.mu)-1
# for snpboostlss, summary for sigma
if(setting$array_id%%4 %in% c(2,3)){
  coef.sigma <- get_coefficients(
    fit$gamma,
    fit$n_step,
    covariates=fit$configs[['covariates']]
  )
  res$coef.sigma <- coef.sigma
  res$sparsity.sigma <- length(coef.sigma)-1
}

# prediction performance
if(setting$array_id%%4 < 2){
  ## MSEP
  pred.snpboost <- predict_snpboost(
    fit,
    new_genotype_file=genotype.pfile,
    new_phenotype_file=phenotype.file,
    phenotype=phenotype,
    subset="test")
  res$MSEP <- pred.snpboost$metric
  ## R2
  res$R2 <- pred.snpboost$Rsquare
  ## loss
  computeLoss <- function (phenotype, prediction.mu, prediction.sigma) {
    rho <- sum(prediction.sigma + (phenotype - prediction.mu)^2 / (2*exp(2*prediction.sigma)))
    return(rho)
  }
  pheno.test <- pred.snpboost$prediction + pred.snpboost$residuals
  sigma2 <- mean(pheno.test^2) - mean(pred.snpboost$prediction^2)
  res$loss  <- computeLoss(
    phenotype = pheno.test,
    prediction.mu = pred.snpboost$prediction,
    prediction.sigma = log(sqrt(sigma2))
  )
  ## PI coverage
  res$PI <- cbind(pred.snpboost$prediction-1.96*sqrt(sigma2),
                  pred.snpboost$prediction+1.96*sqrt(sigma2))
  PI_coverage <- (pheno.test >= res$PI[,1]) & (pheno.test <= res$PI[,2])
  res$PI_cov_rate <- mean(PI_coverage)
  
  # save prediction on test set
  res$pred <- list(mu = pred.snpboost, sigma2 = sigma2)
} else {
  ## MSEP
  pred.snpboostlss <- predict_snpboostlss(
    fit,
    new_genotype_file=genotype.pfile,
    new_phenotype_file=phenotype.file,
    phenotype=phenotype,
    subset="test")
  res$MSEP <- mean((pred.snpboostlss$phenotype - pred.snpboostlss$pred_mu)^2)
  ## R2
  res$R2 <- cor(pred.snpboostlss$phenotype, pred.snpboostlss$pred_mu)^2
  ## loss
  res$loss  <- pred.snpboostlss$metric
  ## PI coverage
  res$PI <- cbind(pred.snpboostlss$pred_mu-1.96*exp(pred.snpboostlss$pred_sigma),
                  pred.snpboostlss$pred_mu+1.96*exp(pred.snpboostlss$pred_sigma))
  PI_coverage <- (pred.snpboostlss$phenotype >= res$PI[,1]) & (pred.snpboostlss$phenotype <= res$PI[,2])
  res$PI_cov_rate <- mean(PI_coverage)
  
  # save prediction on test set
  res$pred <- pred.snpboostlss
}

## computing time
res$time <- run.time

########## Return results ############################

out <- list(setting = setting,
            fit = fit,
            res = res)
save(out,
     file = paste0(dir.res, "/", setting$pheno, "_", setting$method, "_bl.RData"))
