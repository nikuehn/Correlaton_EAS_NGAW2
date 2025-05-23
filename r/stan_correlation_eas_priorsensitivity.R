library(cmdstanr)
library(posterior)
library(bayesplot)
library(matrixStats)
library(tidyverse)


################################################################################
`%notin%` <- Negate(`%in%`)

convert_to_numeric <- function(str) {
  str <- sub("^T", "", str)
  str <- sub("p", ".", str)
  as.numeric(str)
}

inverse_ztransform <- function(rho) {
  return((exp(2 * rho) - 1) / (exp(2 * rho) + 1))
}
ztransform <- function(rho) {
  0.5 * log((1+rho) / (1-rho)) 
}

wid <- 8
asp <- 0.8

################################################################################
dir_data <- file.path('../', 'data/')
dir_stan <- file.path('../', 'stan')
dir_plot <- file.path('../', 'pictures')


################################################################################
# total FAS residuals
totres_fas_combined <- read.csv(file.path(dir_data, 'TotResidAllPer_EAS_ModES_combined.csv'))
dim(totres_fas_combined)

names_target <- names(totres_fas_combined)[!is.na(str_extract(names(totres_fas_combined),pattern = "T[0-9]"))]
freqs_target <- 1/convert_to_numeric(names_target)
n_target <- length(names_target)



################################################################################
# stan models
mod1 <- cmdstan_model(file.path(dir_stan, 'gmm_partition_corrre_cond_miss_tauM_phiM.stan'))
mod_z <- cmdstan_model(file.path(dir_stan, 'gmm_partition_corrre_cond_miss_tauM_phiM_priorz.stan'))



# function to run the stan model
# inputs:
#   k1: index of first target variable
#   k2: index of second target variable
#   model: base for uniform prior on correlation coefficient, otherwise use z-prior
#   pars_z_xx: prior mean and standard deviation for z-transformed correlation coefficient
#   sd_z: if 'given', use value that s provided,
#         otherwise calculate sd for prior from number of records/stations/events
partition_data_rho_eas <- function(k1, k2,
                                   iter_warmup = 500, iter_sampling = 1000, chains = 4,
                                   model = 'base',
                                   pars_z_eq = c(0.55, 0.5),
                                   pars_z_stat = c(0.55, 0.5),
                                   pars_z_rec = c(0.55, 0.5),
                                   sd_z = 'given') {
  data_sel1 <- totres_fas_combined %>%
    select(RSN, EQID, SSN, M, names_target[k1], names_target[k2])
  
  # find target variable which has larger number of data
  n1 <- sum(!is.na(data_sel1[,names_target[k1]]))
  n2 <- sum(!is.na(data_sel1[,names_target[k2]]))
  
  if(n1 > n2) {
    i1 <- 1
    i2 <- 2
    name_1 <- names_target[k1]
    name_2 <- names_target[k2]
  } else {
    i1 <- 2
    i2 <- 1
    name_1 <- names_target[k2]
    name_2 <- names_target[k1]
  }
  data_used <- data_sel1[!is.na(data_sel1[,name_1]),]
  print(dim(data_used))
  
  y_target <- data_used[,c(name_1, name_2)]
  idx_miss_2 <- which(is.na(y_target[,2]))
  idx_obs_2 <- which(!is.na(y_target[,2]))
  
  eq <- as.numeric(factor(data_used$EQID, levels = unique(data_used$EQID)))
  stat <- as.numeric(factor(data_used$SSN, levels = unique(data_used$SSN)))
  
  y_target[idx_miss_2,2] <- -999
  
  data_obs <- data_used[idx_obs_2,]
  if(sd_z != 'given') {
    pars_z_eq[2] <- 2/sqrt(length(unique(data_obs$EQID)) - 3)
    pars_z_stat[2] <- 2/sqrt(length(unique(data_obs$SSN)) - 3)
    pars_z_rec[2] <- 2/sqrt(nrow(data_obs) - 3)
  }
  
  # define predictors for magnitude dependency
  mb <- c(4.5,5.5)
  m1_rec <- if_else(data_used$M <= mb[1], 1,
                    if_else(data_used$M < mb[2], mb[2] - data_used$M, 0))
  m2_rec <- if_else(data_used$M <= mb[1], 0,
                    if_else(data_used$M < mb[2], data_used$M - mb[1], 1))
  
  mageq <- unique(data_used[,c('M','EQID')])$M
  m1_eq <- if_else(mageq <= mb[1], 1,
                   if_else(mageq < mb[2], mb[2] - mageq, 0))
  m2_eq <- if_else(mageq <= mb[1], 0,
                   if_else(mageq < mb[2], mageq - mb[1], 1))
  
  data_list <- list(
    N = nrow(y_target),
    NEQ = max(eq),
    NSTAT = max(stat),
    N_obs_var2 = length(idx_obs_2),
    Y = y_target,
    eq = eq,
    stat = stat,
    idx_obs_var2 = idx_obs_2,
    M1_rec = m1_rec,
    M2_rec = m2_rec,
    M1_eq = m1_eq,
    M2_eq = m2_eq,
    pars_z_eq = pars_z_eq,
    pars_z_stat = pars_z_stat,
    pars_z_rec = pars_z_rec
  )
  
  if(model == 'base') {
    fit_partition <- mod1$sample(
      data = data_list,
      seed = 1701,
      parallel_chains = 2,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = 50,
      show_exceptions = FALSE
    )
  } else {
    fit_partition <- mod_z$sample(
      data = data_list,
      seed = 1701,
      parallel_chains = 2,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = 50,
      show_exceptions = FALSE
    )
  }
  print(fit_partition$cmdstan_diagnose())
  print(fit_partition$diagnostic_summary())
  draws <- fit_partition$draws()
  rv <- as_draws_rvars(draws)
  
  rv$rho_total_sm <- (rv$phi_ss_sm[1] * rv$phi_ss_sm[2] * rv$rho_rec +
                        rv$phi_s2s[1] * rv$phi_s2s[2] * rv$rho_stat +
                        rv$tau_sm[1] * rv$tau_sm[2] * rv$rho_eq) /
    (sqrt(rv$phi_ss_sm[1]^2 + rv$phi_s2s[1]^2 + rv$tau_sm[1]^2) *
       sqrt(rv$phi_ss_sm[2]^2 + rv$phi_s2s[2]^2 + rv$tau_sm[2]^2))
  
  rv$rho_total_lm <- (rv$phi_ss_lm[1] * rv$phi_ss_lm[2] * rv$rho_rec +
                        rv$phi_s2s[1] * rv$phi_s2s[2] * rv$rho_stat +
                        rv$tau_lm[1] * rv$tau_lm[2] * rv$rho_eq) /
    (sqrt(rv$phi_ss_lm[1]^2 + rv$phi_s2s[1]^2 + rv$tau_lm[1]^2) *
       sqrt(rv$phi_ss_lm[2]^2 + rv$phi_s2s[2]^2 + rv$tau_lm[2]^2))
  
  tmp <- data.frame(target1 = name_1,
                    target2 = name_2,
                    rho_eq = rv$rho_eq,
                    rho_stat = rv$rho_stat,
                    rho_rec = rv$rho_rec,
                    rho_total_sm = rv$rho_total_sm, rho_total_lm = rv$rho_total_lm,
                    tau_target1_sm = rv$tau_sm[1], tau_target2_sm = rv$tau_sm[2],
                    tau_target1_lm = rv$tau_lm[1], tau_target2_lm = rv$tau_lm[2],
                    phi_s2s_target1 = rv$phi_s2s[1], phi_s2s_target2 = rv$phi_s2s[2],
                    phi_ss_target1_sm = rv$phi_ss_sm[1], phi_ss_target2_sm = rv$phi_ss_sm[2],
                    phi_ss_target1_lm = rv$phi_ss_lm[1], phi_ss_target2_lm = rv$phi_ss_lm[2])
  
  tmp2 <- summarise_draws(
    subset(draws, variable = c('phi','tau','rho','z_'), regex = TRUE),
    "mean",
    "median",
    "sd",
    "mad",
    "rhat",
    "ess_bulk",
    "ess_tail",
    ~quantile(.x, probs = c(0.025, 0.975))
  ) %>%
    set_names("variable","mean","median","sd","mad","rhat","ess_bulk","ess_tail","q2.5","q97.5")
  
  return(list(rv = tmp, summary = tmp2, data = data_list))
}


# run stan models for frequencies 1 and 5Hz
# for time, run only 20 warmup and 200 post warmup
k1 <- 14 # 1Hz
k2 <- 8
print(c(names_target[k1], names_target[k2]))


n_sampling <- 200
n_warmup <- 200
n_chains <- 4
n_post <- n_chains * n_sampling
result_base <- partition_data_rho_eas(k1, k2, iter_warmup = n_warmup,
                                      iter_sampling = n_sampling, chains = n_chains, model = 'base')
result_z <- partition_data_rho_eas(k1, k2, iter_warmup = n_warmup,
                                   iter_sampling = n_sampling, chains = n_chains, model = 'z')

# informative prior from BA19
cor_ba19 <- read.csv(file.path(dir_data, "cor_ba19.csv"))

z_eq <- ztransform(cor_ba19$eq[k2])
z_stat <-  ztransform(cor_ba19$stat[k2])
z_rec <-  ztransform(cor_ba19$rec[k2])

result_z_inf <- partition_data_rho_eas(k1, k2, iter_warmup = n_warmup,
                                       iter_sampling = n_sampling, chains = n_chains, model = 'z',
                                       pars_z_eq = c(z_eq,NA), pars_z_stat = c(z_stat,NA), pars_z_rec = c(z_rec,NA),
                                       sd_z = 'calc')


(pl <- patchwork::wrap_plots(
  data.frame(unif_posterior = as_draws_matrix(result_base$rv$rho_eq),
             prior1_posterior = as_draws_matrix(result_z$rv$rho_eq),
             prior2_posterior = as_draws_matrix(result_z_inf$rv$rho_eq),
             prior1_prior = tanh(rnorm(n_post, mean = 0.55, sd  = 0.5)),
             prior2_prior = tanh(rnorm(n_post,
                                       mean = result_z_inf$data$pars_z_eq[1],
                                       sd  = result_z_inf$data$pars_z_eq[2]))
  ) %>% set_names(c('unif_posterior','prior1_posterior','prior2_posterior','prior1_prior','prior2_prior')) %>%
    pivot_longer(everything(), names_to = c('model','type'), names_sep = '_') %>%
    ggplot(aes(x = value, color = model, linetype = type)) +
    geom_density(linewidth = 1.5) +
    scale_color_manual(values = c('red','blue', 'black')) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', linewidth = 1.5) +
    labs(x = 'rho_dB')
  ,
  data.frame(unif_posterior = as_draws_matrix(result_base$rv$rho_stat),
             prior1_posterior = as_draws_matrix(result_z$rv$rho_stat),
             prior2_posterior = as_draws_matrix(result_z_inf$rv$rho_stat),
             prior1_prior = tanh(rnorm(n_post, mean = 0.55, sd  = 0.5)),
             prior2_prior = tanh(rnorm(n_post,
                                       mean = result_z_inf$data$pars_z_stat[1],
                                       sd  = result_z_inf$data$pars_z_stat[2]))
  ) %>% set_names(c('unif_posterior','prior1_posterior','prior2_posterior','prior1_prior','prior2_prior')) %>%
    pivot_longer(everything(), names_to = c('model','type'), names_sep = '_') %>%
    ggplot(aes(x = value, color = model, linetype = type)) +
    geom_density(linewidth = 1.5) +
    scale_color_manual(values = c('red','blue', 'black')) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', linewidth = 1.5) +
    labs(x = 'rho_dS')
  ,
  data.frame(unif_posterior = as_draws_matrix(result_base$rv$rho_rec),
             prior1_posterior = as_draws_matrix(result_z$rv$rho_rec),
             prior2_posterior = as_draws_matrix(result_z_inf$rv$rho_rec),
             prior1_prior = tanh(rnorm(n_post, mean = 0.55, sd  = 0.5)),
             prior2_prior = tanh(rnorm(n_post,
                                       mean = result_z_inf$data$pars_z_rec[1],
                                       sd  = result_z_inf$data$pars_z_rec[2]))
  ) %>% set_names(c('unif_posterior','prior1_posterior','prior2_posterior','prior1_prior','prior2_prior')) %>%
    pivot_longer(everything(), names_to = c('model','type'), names_sep = '_') %>%
    ggplot(aes(x = value, color = model, linetype = type)) +
    geom_density(linewidth = 1.5) +
    scale_color_manual(values = c('red','blue', 'black')) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', linewidth = 1.5) +
    labs(x = 'rho_dWS')
  , ncol = 2) + 
    patchwork::plot_annotation(title = 'Correlation between 1Hz and 5Hz')
)

ggsave(file.path(dir_plot, sprintf('plot_priorsensitivity.png')), pl,
       width = 2 * wid, height = 2 * asp * wid)
