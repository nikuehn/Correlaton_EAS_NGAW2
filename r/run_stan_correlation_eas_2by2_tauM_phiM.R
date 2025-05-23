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
# define Stan model to partition residuals and calculate correlations
mod <- cmdstan_model(file.path(dir_stan, 'gmm_partition_corrre_cond_miss_tauM_phiM.stan'))

# function to run the stan model
# inputs:
#   k1: index of first target variable
#   k2: index of second target variable
partition_data_rho_eas <- function(k1, k2,
                                   iter_warmup = 500, iter_sampling = 1000, chains = 4) {
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
    M2_eq = m2_eq
  )
  
  fit_partition <- mod$sample(
    data = data_list,
    seed = 1701,
    parallel_chains = 2,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 50,
    show_exceptions = FALSE
  )
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
    subset(draws, variable = c('phi','tau','rho'), regex = TRUE),
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
  
  return(list(rv = tmp, summary = tmp2))
}


# run stan model for frequencies 1 and 5Hz
k1 <- 8
k2 <- 14
print(freqs_target[c(k1,k2)])
result_stan <- partition_data_rho_eas(k1,k2)

print(result_stan$summary)

# plot posteroir distributions of correlation coefficients
(pl <- data.frame(type = c('rho_dB','rho_dS','rho_dWS'),
           rho = c(result_stan$rv$rho_eq, result_stan$rv$rho_stat, result_stan$rv$rho_rec)) |>
  ggplot(aes(y=type, xdist = rho)) +
  ggdist::stat_slabinterval() +
  labs(title = paste0('Correlation between F = ',freqs_target[k1],'Hz and F = ',freqs_target[k2],'Hz'))
)
ggsave(file.path(dir_plot, sprintf('plot_rho_5Hz1Hz.png')), pl,
       width = wid, height = asp * wid)





# run stan model for all frequency combinations 
results_cor_rv <- list()
results_cor_summary <- list()
k <- 1
for(k1 in 1:n_target) {
  for(k2 in (k1+1):n_target) {
    if(k2 > length(names_target)) {break}
    
    print(c(names_target[k1], names_target[k2]))
    tmp <- partition_data_rho_eas(k1, k2, iter_warmup = 500, iter_sampling = 1000, chains = 4)
    results_cor_rv[[k]] <- tmp$rv
    results_cor_summary[[k]] <- tmp$summary
    k <- k + 1
  }
}


