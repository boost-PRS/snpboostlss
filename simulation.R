library(data.table)
library(parallel)
library(bigsnpr)
library(ggplot2)
library(cowplot)
library(readr)
library(gamlss)
library(pgenlibr)
library(IDSL.IPA)
library(forstringr)
library(latex2exp)
library(tidyverse)

settings <- expand_grid(
  n = 20000,
  s = c(0.001, 0.01),
  h2 = c(0.1, 0.7)
)

array.id <- 1 # array.id represent the simulation number

setting.idx <- array.id %/% 100 + 1
sim <- array.id %% 100

n <- settings$n[setting.idx]
p <- n
s <- settings$s[setting.idx]
h2 <- settings$h2[setting.idx]

set.seed(2024*sim + setting.idx * 10)


dir.sim <- "../simulation"
bim <- fread(paste0(dir.sim, "/synthetic_v1_chr-22_filtered.bim"))  # HAPNEST data chromosome 22
fam <- fread(paste0(dir.sim, "/synthetic_v1_chr-22_filtered.fam"))  # HAPNEST data chromosome 22
split_default <- c(rep("train",0.5*n),
                   rep("val",0.2*n), 
                   rep("test",0.3*n))

# absolute path of working directory for current simulation run
dir <- paste0(dir.sim, "/n",n,"_s",s,"_heri",h2,"_sim_", sim)
if (file.exists(dir)){
  setwd(dir)
} else {
  dir.create(dir)
  setwd(dir)
}

# set directory to store results
dir.res <- paste0(dir.sim, "/results")
if (!file.exists(dir.res)){
  dir.create(dir.res)
} 

########## Generate genotype #########################
DataGen.start.time <- Sys.time()

# generate n observations
people <- sample(1:nrow(fam),n)
people <- sort(people)
sample_fam <- fam[people,]
# write a text file containing the names of selected n observations
write.table(sample_fam[,1:2], file = paste0(dir,"/people_n",n,"_s",s,"_heri",h2,"_sim_", sim, ".txt"), sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# split n observations into 50% training, 20% validation and 30% test sets
pheno_simulated <- data.frame(FID= sample_fam$V1, 
                              IID = sample_fam$V2, 
                              split=sample(split_default, length(split_default), replace = FALSE))

# generate p variants
variants <- sample(1:nrow(bim),p)
variants <- sort(variants)
sample_bim <- bim[variants,]
# write a text file containing the names of selected p variants
write.table(sample_bim[,2], file = paste0(dir,"/variants_n",n,"_s",s,"_heri",h2,"_sim_", sim,".txt"), sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# get the subsetted bfiles containing p variants and n observations
system(paste0(dir.sim, "/plink --bfile ", dir.sim, "/synthetic_v1_chr-22_filtered --keep ",
              dir, "/people_n", n,"_s",s,"_heri",h2,"_sim_", sim, 
              ".txt --extract ", dir, "/variants_n", n,"_s",s,"_heri",h2,"_sim_",sim, 
              ".txt --make-bed --out geno_n",n,"_s",s,"_heri",h2,"_sim_", sim))

# transform bfiles (file format for PLINK1) to pfiles (file formate for PLINK2)
system(paste0(dir.sim, "/plink2 --bfile geno_n",n,"_s",s,"_heri",h2,"_sim_",sim,
              " --make-pgen vzs --out geno_n",n,"_s",s,"_heri",h2,"_sim_",sim))
system(paste0(dir.sim, "/plink2 --bfile geno_n",n,"_s",s,"_heri",h2,"_sim_",sim,
              " --make-pgen --out geno_n",n,"_s",s,"_heri",h2,"_sim_",sim))


########## Generate gamma ############################
# select informative variants for sigma
variants_infor_index_sigma <- sample(1:nrow(sample_bim),p*s)
variants_infor_index_sigma <- sort(variants_infor_index_sigma)
variants_infor_sigma <- data.frame(rsID = sample_bim$V2[variants_infor_index_sigma], 
                                   CHR = sample_bim$V1[variants_infor_index_sigma], 
                                   Pos = sample_bim$V4[variants_infor_index_sigma],
                                   ALT = sample_bim$V5[variants_infor_index_sigma])
# generate gamma from U(-0.25, 0.25)
variants_infor_sigma$gamma <- runif(length(variants_infor_index_sigma), min = -0.25, max = 0.25)
# save the true coefficients 
fwrite(variants_infor_sigma,
       paste0(dir,"/true_gamma_n",n,"_s",s,"_heri",h2,"_sim_", sim,".txt"), 
       sep="\t")

# calculate linear predictor for scale (log_sigma) using PLINK command
system(paste0(dir.sim, "/plink2 --pfile geno_n",n,"_s",s,"_heri",h2,"_sim_",sim,
              " --score true_gamma_n",n,"_s",s,"_heri",h2,"_sim_",sim,
              ".txt 1 4 5 cols=scoresums --out log_sigma_n",n,"_s",s,"_heri",h2,"_sim_", sim
))



########## Generate beta #############################
# select informative variants for mu
variants_infor_index_mu <- sample(1:nrow(sample_bim),p*s)
variants_infor_index_mu <- sort(variants_infor_index_mu)
variants_infor_mu <- data.frame(rsID = sample_bim$V2[variants_infor_index_mu], 
                                CHR = sample_bim$V1[variants_infor_index_mu], 
                                Pos = sample_bim$V4[variants_infor_index_mu],
                                ALT = sample_bim$V5[variants_infor_index_mu])
# read in calculated linear predictor for scale
log_sigma <- fread(paste0(dir,"/log_sigma_n",n,"_s",s,"_heri",h2,"_sim_",sim,".sscore"))
log_sigma <- log_sigma %>% rename(IID = `#IID`, log_sigma_true = SCORE1_SUM) 
pheno_simulated <- left_join(pheno_simulated, log_sigma, by = "IID")
sigma <- exp(log_sigma$log_sigma_true)
# calculate average error term variance
avg_sigma2 <- mean(sigma^2)
# generate beta from normal distribution based on average error term variance
variants_infor_mu$beta <- rnorm(length(variants_infor_index_mu), mean = 0, sd = sqrt((avg_sigma2*h2/(1-h2))/length(variants_infor_index_mu)))
# save the true coefficients 
fwrite(variants_infor_mu,paste0(dir, "/true_beta_n",n,"_s",s,"_heri",h2,"_sim_",sim,".txt"), sep="\t")

## calculate linear predictor for location (eta_mu)
# load pgen data
genotype.pfile <- paste0("geno_n",n,"_s",s,"_heri",h2,"_sim_",sim)
pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar'))
pgen.tmp <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar)
pgenlibr::ClosePvar(pvar)

# subset pgen data for informative variants for mu
buf <- pgenlibr::ReadList(pgen.tmp, variants_infor_index_mu, meanimpute=F)
variants.mu <- as.data.table(buf)

# calculate mu as product of standardized genotype data and beta
mu <- scale(as.matrix(variants.mu)) %*% variants_infor_mu$beta
# format mu into data frame with IID
iids.psam <- data.table::fread(paste0(genotype.pfile, '.psam'))
mu <- data.frame(IID = iids.psam$IID, mu_true = mu[,1])

########## Generate phenotype ########################
# read in calculated linear predictor for mu
pheno_simulated <- left_join(pheno_simulated, mu, by = "IID")
# generate phenotype from normal distribution
pheno_simulated$pheno=rnorm(nrow(sample_fam), mean = mu$mu_true, sd = sigma) 
# calculate h2 for generated data 
h2_actual <- var(pheno_simulated$mu_true)/var(pheno_simulated$pheno)
# calculate total variance and residual variance
var <- pheno_simulated %>%
  filter(split == "test") %>%
  mutate(error = pheno - mu_true) %>%
  summarise(
    var.y = var(pheno),
    var.e = var(error)
  )
var.y <- var$var.y
var.e <- var$var.e

fwrite(pheno_simulated,paste0(dir,"/pheno_n",n,"_s",s,"_heri",h2,"_sim_",sim,".phe"), sep=" ")

remove(people, variants,
       variants_infor_index_sigma,
       variants_infor_index_mu, log_sigma,
       sigma,avg_sigma2,mu)

remove(bim, buf, fam, iids.psam, pgen.tmp, pvar, variants.mu)

# file.remove(paste0("geno_n",n,"_s",s,"_heri",h2,"_sim_", sim, ".pvar"))

DataGen.end.time <- Sys.time()

DataGen.time <- difftime(DataGen.end.time, DataGen.start.time, unit="mins")

paste("Data generation DONE. It takes", round(DataGen.time,4), units(DataGen.time))



########## Fit snpboost FSL ##############################
### load functions
source("../snpboost/lm_boost.R")
source("../snpboost/snpboost.R")
source("../snpboost/functions_snpboost.R")
source("../snpboost/functions_snpnet.R")
# get the snpboost functions from https://github.com/boost-PRS/snpboost


genotype.pfile <- paste0("geno_n",n,"_s",s,"_heri",h2,"_sim_",sim)  # genotype file
phenotype.file <- paste0("pheno_n",n,"_s",s,"_heri",h2,"_sim_",sim,".phe") # phenotype file, .phe, tabular data
phenotype <- "pheno" # name of phenotype column in phenotype file
covariates <- NULL # name of covariates for mu

configs <- list(plink2.path = paste0(dir.sim,"/plink2"), # path to plink file
                zstdcat.path = "zstd", # path to zstdcat (or zstd)
                zcat.path = "cat",
                results.dir = paste0(dir,"/results"), # results folder
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


snpboostFSL.start.time <- Sys.time()

fit_snpboost_fsl <- snpboost(
  genotype.pfile = genotype.pfile,
  phenotype.file = phenotype.file,
  phenotype = phenotype,
  covariates= covariates,
  configs = configs,
  split.col = "split", ### name of column in phenotype file indication train, val, test set
  p_batch = 1000, # batch size
  m_batch = 1000, # maximum number of boosting steps per batch
  b_max = 20000, # maximum number of batches
  b_stop = 2, # stop if performance not improving for more than 2 batches
  sl= 0.1, # learning rate (default=0.1)
  give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
  family = "gaussian", # family argument: currently only gaussian
  metric = "MSEP" # loss function/metric: currently only negative log-likelihood for gaussian
)

snpboostFSL.end.time <- Sys.time()

snpboostFSL.time <- difftime(snpboostFSL.end.time, snpboostFSL.start.time, unit="mins")

paste("snpboostFSL DONE. It takes", round(snpboostFSL.time,4), units(snpboostFSL.time))

# calculate performance measures
res_snpboost_fsl <- list()
## number of included variants
coef_snpboost_fsl.mu <- get_coefficients(fit_snpboost_fsl$beta,
                                     fit_snpboost_fsl$n_step,
                                     covariates=fit_snpboost_fsl$configs[['covariates']])
res_snpboost_fsl$coef <- coef_snpboost_fsl.mu
res_snpboost_fsl$sparsity <- length(coef_snpboost_fsl.mu)-1
## true positive rate
TP_snpboost_fsl <- intersect(variants_infor_mu$rsID, substring(names(coef_snpboost_fsl.mu),1,18))
res_snpboost_fsl$TPR <- length(TP_snpboost_fsl)/nrow(variants_infor_mu)
## precision
res_snpboost_fsl$precision <- length(TP_snpboost_fsl)/length(coef_snpboost_fsl.mu)
## MSEP
pred.snpboost_fsl <- predict_snpboost(
  fit_snpboost_fsl,
  new_genotype_file=genotype.pfile,
  new_phenotype_file=phenotype.file,
  phenotype=phenotype,
  subset="test")
res_snpboost_fsl$MSEP <- pred.snpboost_fsl$metric
## R2
res_snpboost_fsl$R2 <- pred.snpboost_fsl$Rsquare
## loss
computeLoss <- function (phenotype, prediction.mu, prediction.sigma) {
  rho <- sum(prediction.sigma + (phenotype - prediction.mu)^2 / (2*exp(2*prediction.sigma)))
  return(rho)
}
sigma2_fsl <- sum(pred.snpboost_fsl$residuals^2)/(length(pred.snpboost_fsl$prediction)-length(coef_snpboost_fsl.mu))
res_snpboost_fsl$loss  <- computeLoss(
  phenotype = pred.snpboost_fsl$prediction + pred.snpboost_fsl$residuals,
  prediction.mu = pred.snpboost_fsl$prediction,
  prediction.sigma = log(sqrt(sigma2_fsl)))
## PI coverage
res_snpboost_fsl$PI <- cbind(pred.snpboost_fsl$prediction-1.96*sqrt(sigma2_fsl),
                         pred.snpboost_fsl$prediction+1.96*sqrt(sigma2_fsl))
PI_coverage_snpboost_fsl <- ((pred.snpboost_fsl$prediction+pred.snpboost_fsl$residuals) >= res_snpboost_fsl$PI[,1]) &
               ((pred.snpboost_fsl$prediction+pred.snpboost_fsl$residuals) <= res_snpboost_fsl$PI[,2])
res_snpboost_fsl$PI_cov_rate <- mean(PI_coverage_snpboost_fsl)
## computing time
res_snpboost_fsl$time <- snpboostFSL.time


########## Fit snpboostlss with fixed sl ##################

## load functions
source("../snpboostlss/lmlss_boost.R")
source("../snpboostlss/snpboostlss.R")
source("../snpboostlss/functions_snpboostlss.R")
source("../snpboostlss/functions_snpboost.R")
source("../snpboostlss/functions_snpnet.R")
# get the snpboost functions from https://github.com/boost-PRS/snpboostlss

genotype.pfile <- paste0("geno_n",n,"_s",s,"_heri",h2,"_sim_",sim)  # genotype file
phenotype.file <- paste0("pheno_n",n,"_s",s,"_heri",h2,"_sim_",sim,".phe") # phenotype file, .phe, tabular data
phenotype <- "pheno" # name of phenotype column in phenotype file
covariates.mu <- NULL # name of covariates for mu
covariates.sigma <- NULL # name of covariates for sigma

configs <- list(plink2.path = paste0(dir.sim,"/plink2"), # path to plink file
                zcat.path = "cat",
                # sed.path = "../GnuWin32/bin/sed", # for running on Windows
                results.dir = paste0(dir,"/results"),
                save = FALSE, ## if true, results are saved after each batch
                prevIter = 0,
                missing.rate = 0.1,
                MAF.thresh = 0.001,
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

snpboostlssFSL.start.time <- Sys.time()

fit_snpboostlss_fsl <- snpboostlss(
  data.dir = paste0(dir,"/"),
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
  sl= 0.1, # learning rate (default=NULL, i.e. adaptive sl).
  # A value different from NULL means using the same fixed sl for mu and sigma.
  give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
  family = "gaussian", # family argument: currently only gaussian
  metric = "loss" # loss function/metric: currently only negative log-likelihood for gaussian
)

snpboostlssFSL.end.time <- Sys.time()

snpboostlssFSL.time <- difftime(snpboostlssFSL.end.time, snpboostlssFSL.start.time, unit="mins")

paste("snpboostlss with fixed step length DONE. It takes", snpboostlssFSL.time, units(snpboostlssFSL.time))

# calculate performance measures
res_snpboostlss_fsl <- list()
## number of included variants
coef_snpboostlss_fsl.mu <- get_coefficients(
  fit_snpboostlss_fsl$beta,
  fit_snpboostlss_fsl$n_step,
  covariates=fit_snpboostlss_fsl$configs[['covariates']])
res_snpboostlss_fsl$coef.mu <- coef_snpboostlss_fsl.mu
res_snpboostlss_fsl$sparsity.mu <- length(coef_snpboostlss_fsl.mu)-1
coef_snpboostlss_fsl.sigma <- get_coefficients(
  fit_snpboostlss_fsl$gamma,
  fit_snpboostlss_fsl$n_step,
  covariates=fit_snpboostlss_fsl$configs[['covariates']])
res_snpboostlss_fsl$coef.sigma <- coef_snpboostlss_fsl.sigma
res_snpboostlss_fsl$sparsity.sigma <- length(coef_snpboostlss_fsl.sigma)-1
## true positive rate
TP_snpboostlss_fsl.mu <- intersect(variants_infor_mu$rsID, substring(names(coef_snpboostlss_fsl.mu),1,18))
res_snpboostlss_fsl$TPR.mu <- length(TP_snpboostlss_fsl.mu)/nrow(variants_infor_mu)
TP_snpboostlss_fsl.sigma <- intersect(variants_infor_sigma$rsID, substring(names(coef_snpboostlss_fsl.sigma),1,18))
res_snpboostlss_fsl$TPR.sigma <- length(TP_snpboostlss_fsl.sigma)/nrow(variants_infor_sigma)
## precision
res_snpboostlss_fsl$precision.mu <- length(TP_snpboostlss_fsl.mu)/length(coef_snpboostlss_fsl.mu)
res_snpboostlss_fsl$precision.sigma <- length(TP_snpboostlss_fsl.sigma)/length(coef_snpboostlss_fsl.sigma)
## MSEP
pred.snpboostlss.fsl <- predict_snpboostlss(
  fit_snpboostlss_fsl,
  new_genotype_file=genotype.pfile,
  new_phenotype_file=phenotype.file,
  phenotype=phenotype,
  subset="test")
res_snpboostlss_fsl$MSEP <- mean((pred.snpboostlss.fsl$phenotype - pred.snpboostlss.fsl$pred_mu)^2)
## R2
res_snpboostlss_fsl$R2 <- cor(pred.snpboostlss.fsl$phenotype, pred.snpboostlss.fsl$pred_mu)^2
## loss
res_snpboostlss_fsl$loss  <- pred.snpboostlss.fsl$metric
## PI coverage
res_snpboostlss_fsl$PI <- cbind(pred.snpboostlss.fsl$pred_mu-1.96*exp(pred.snpboostlss.fsl$pred_sigma),
                                pred.snpboostlss.fsl$pred_mu+1.96*exp(pred.snpboostlss.fsl$pred_sigma))
PI_coverage_snpboostlss_fsl <- (pred.snpboostlss.fsl$phenotype >= res_snpboostlss_fsl$PI[,1]) &
               (pred.snpboostlss.fsl$phenotype <= res_snpboostlss_fsl$PI[,2])
res_snpboostlss_fsl$PI_cov_rate <- mean(PI_coverage_snpboostlss_fsl)
## computing time
res_snpboostlss_fsl$time <- snpboostlssFSL.time

########## Fit snpboostlss with adaptive sl ##################

source("/home/wuq/RESEARCH/project1/algorithm/lmlss_boost.R")
source("/home/wuq/RESEARCH/project1/algorithm/snpboostlss.R")
source("/home/wuq/RESEARCH/project1/algorithm/functions_snpboostlss.R")
source("/home/wuq/RESEARCH/project1/algorithm/functions_snpboost.R")
source("/home/wuq/RESEARCH/project1/algorithm/functions_snpnet.R")

genotype.pfile <- paste0("geno_n",n,"_s",s,"_heri",h2,"_sim_",sim)  # genotype file
phenotype.file <- paste0("pheno_n",n,"_s",s,"_heri",h2,"_sim_",sim,".phe") # phenotype file, .phe, tabular data
phenotype <- "pheno" # name of phenotype column in phenotype file
covariates.mu <- NULL # name of covariates for mu
covariates.sigma <- NULL # name of covariates for sigma

configs <- list(plink2.path = paste0(dir.sim,"/plink2"), # path to plink file
                zcat.path = "cat",
                # sed.path = "../GnuWin32/bin/sed", # for running on Windows
                results.dir = paste0(dir,"/results"),
                save = FALSE, ## if true, results are saved after each batch
                prevIter = 0,
                missing.rate = 0.1,
                MAF.thresh = 0.001,
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

snpboostlssASL.start.time <- Sys.time()

fit_snpboostlss_asl <- snpboostlss(
  data.dir = paste0(dir,"/"),
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
  sl= NULL, # learning rate (default=NULL, i.e. adaptive sl).
  # A value different from NULL means using the same fixed sl for mu and sigma.
  give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
  family = "gaussian", # family argument: currently only gaussian
  metric = "loss" # loss function/metric: currently only negative log-likelihood for gaussian
)

snpboostlssASL.end.time <- Sys.time()

snpboostlssASL.time <- difftime(snpboostlssASL.end.time, snpboostlssASL.start.time, unit="mins")

paste("snpboostlss with adaptive step length DONE. It takes", snpboostlssASL.time, units(snpboostlssASL.time))

# calculate performance measures
res_snpboostlss_asl <- list()
## number of included variants
coef_snpboostlss_asl.mu <- get_coefficients(
  fit_snpboostlss_asl$beta,
  fit_snpboostlss_asl$n_step,
  covariates=fit_snpboostlss_asl$configs[['covariates']])
res_snpboostlss_asl$coef.mu <- coef_snpboostlss_asl.mu
res_snpboostlss_asl$sparsity.mu <- length(coef_snpboostlss_asl.mu)-1
coef_snpboostlss_asl.sigma <- get_coefficients(
  fit_snpboostlss_asl$gamma,
  fit_snpboostlss_asl$n_step,
  covariates=fit_snpboostlss_asl$configs[['covariates']])
res_snpboostlss_asl$coef.sigma <- coef_snpboostlss_asl.sigma
res_snpboostlss_asl$sparsity.sigma <- length(coef_snpboostlss_asl.sigma)-1
## true positive rate
TP_snpboostlss_asl.mu <- intersect(variants_infor_mu$rsID, substring(names(coef_snpboostlss_asl.mu),1,18))
res_snpboostlss_asl$TPR.mu <- length(TP_snpboostlss_asl.mu)/nrow(variants_infor_mu)
TP_snpboostlss_asl.sigma <- intersect(variants_infor_sigma$rsID, substring(names(coef_snpboostlss_asl.sigma),1,18))
res_snpboostlss_asl$TPR.sigma <- length(TP_snpboostlss_asl.sigma)/nrow(variants_infor_sigma)
## precision
res_snpboostlss_asl$precision.mu <- length(TP_snpboostlss_asl.mu)/length(coef_snpboostlss_asl.mu)
res_snpboostlss_asl$precision.sigma <- length(TP_snpboostlss_asl.sigma)/length(coef_snpboostlss_asl.sigma)
## MSEP
pred.snpboostlss.asl <- predict_snpboostlss(
  fit_snpboostlss_asl,
  new_genotype_file=genotype.pfile,
  new_phenotype_file=phenotype.file,
  phenotype=phenotype,
  subset="test")
res_snpboostlss_asl$MSEP <- mean((pred.snpboostlss.asl$phenotype - pred.snpboostlss.asl$pred_mu)^2)
## R2
res_snpboostlss_asl$R2 <- cor(pred.snpboostlss.asl$phenotype, pred.snpboostlss.asl$pred_mu)^2
## loss
res_snpboostlss_asl$loss  <- pred.snpboostlss.asl$metric
## PI coverage
res_snpboostlss_asl$PI <- cbind(pred.snpboostlss.asl$pred_mu-1.96*exp(pred.snpboostlss.asl$pred_sigma),
                                pred.snpboostlss.asl$pred_mu+1.96*exp(pred.snpboostlss.asl$pred_sigma))
PI_coverage_snpboostlss_asl <- (pred.snpboostlss.asl$phenotype >= res_snpboostlss_asl$PI[,1]) &
  (pred.snpboostlss.asl$phenotype <= res_snpboostlss_asl$PI[,2])
res_snpboostlss_asl$PI_cov_rate <- mean(PI_coverage_snpboostlss_asl)
## computing time
res_snpboostlss_asl$time <- snpboostlssASL.time

########## Return results ############################

out <- list(sample_fam = sample_fam,
            sample_bim = sample_bim,
            variants_infor_mu = variants_infor_mu,
            variants_infor_sigma = variants_infor_sigma,
            pheno_simulated = pheno_simulated,
            h2_actual = h2_actual,
            var.y = var.y,
            var.e = var.e,
            DataGen.time = DataGen.time,
            fit_snpboost_fsl = fit_snpboost_fsl,
            res_snpboost_fsl = res_snpboost_fsl,
            fit_snpboostlss_fsl = fit_snpboostlss_fsl,
            res_snpboostlss_fsl = res_snpboostlss_fsl,
            fit_snpboostlss_asl = fit_snpboostlss_asl,
            res_snpboostlss_asl = res_snpboostlss_asl
)

assign(paste0("res_n",n,"_s",s,"_heri",h2,"_sim_",sim), out)
save(list = c(paste0("res_n",n,"_s",s,"_heri",h2,"_sim_",sim)),
     file = paste0(dir.res, "/results_n",n,"_s",s,"_heri",h2,"_sim_",sim,".RData"))



##########################################################################################################

######################### Longitudinal data generation #########################

set.seed(1234)

# max number of longitudinal points
t.max <- 100 

s <- s
heri <- h2

# create a data frame to store information about longitudinal data
longi <- data.frame(matrix(NA, ncol = (2*t.max+8), nrow = 6000))
colnames(longi) <- c("id","s","heri","sim","mu_true","sigma_true",
                     "PI_low","PI_up","sigma_hat_bl",
                     paste0("y_",0:(t.max-1)),
                     paste0("sigma_hat_longi_",2:t.max))

# fill in some basic information
longi$id <- out$pheno_simulated$IID[out$pheno_simulated$split=="test"]
longi$s <- s
longi$heri <- heri
longi$sim <- sim
longi$mu_true <- out$pheno_simulated$mu_true[out$pheno_simulated$split=="test"]
longi$sigma_true <- exp(out$pheno_simulated$log_sigma_true[out$pheno_simulated$split=="test"])
longi$PI_low <- out$res_snpboostlss_asl$PI[,1]
longi$PI_up <- out$res_snpboostlss_asl$PI[,2]
longi$sigma_hat_bl <- (longi$PI_up - longi$PI_low)/2/1.96
longi$y_0 <- out$pheno_simulated$pheno[out$pheno_simulated$split=="test"]

for(i in 1:6000){
  # generate longitudinal data
  longi[i, paste0("y_", 1:(t.max-1))] <- rnorm(n=(t.max-1),
                                               mean=longi$mu_true[i],
                                               sd=longi$sigma_true[i])
  # calculate within person variability
  for(j in 2:t.max){
    longi[i, paste0("sigma_hat_longi_",j)] <- sd(longi[i,paste0("y_",0:(j-1))])
  }
}


# calculate prediction performance for sigma estimators
corr_sigma <- apply(longi[,c("sigma_hat_bl", paste0("sigma_hat_longi_", 2:t.max))],
                  2, # by column
                  function(hat){cor(hat, longi$sigma_true)})
# save results
out$longi <- longi
out$corr_sigma <- corr_sigma
assign(paste0("res_n",n,"_s",s,"_heri",h2,"_sim_",sim), out)
save(list = c(paste0("res_n",n,"_s",s,"_heri",h2,"_sim_",sim)),
     file = paste0(dir.res, "/results_n",n,"_s",s,"_heri",h2,"_sim_",sim,".RData"))



