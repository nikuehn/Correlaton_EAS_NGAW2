/* ***********************************************
 Model that partitions total resduals of multiple target variables
 into even terms, site terms, and singe-site-residuals, 
 and estmates correlations between them.
 
 This model is for two target variables, and the second is
 modeled contional on the first.

 This model assumes that the first target variable is complete,
 and the second variable can be incomplete.

 *********************************************** */

data {
  int<lower=1> N;
  int<lower=1> NEQ;
  int<lower=1> NSTAT;
  int<lower=1,upper=N> N_obs_var2;

  matrix[N, 2] Y;    // log psa values - 1 is PGA

  array[N] int<lower=1,upper=NEQ> eq;
  array[N] int<lower=1,upper=NSTAT> stat;

  array[N_obs_var2] int<lower=1,upper=N> idx_obs_var2;

  vector[NEQ] M1_eq;
  vector[NEQ] M2_eq;
  vector[N] M1_rec;
  vector[N] M2_rec;

  vector[2] pars_z_eq;
  vector[2] pars_z_stat;
  vector[2] pars_z_rec;
}

transformed data {
  int N_target = 2;
  vector[N_target] zerovec = rep_vector(0,N_target);
}

parameters {
  vector[N_target] c0;

  vector<lower=0>[N_target] phi_ss_sm;
  vector<lower=0>[N_target] phi_ss_lm;
  vector<lower=0>[N_target] phi_s2s;
  vector<lower=0>[N_target] tau_sm;
  vector<lower=0>[N_target] tau_lm;

  real z_rec;
  real z_eq;
  real z_stat;

  matrix[NEQ, N_target] deltaB;
  matrix[NSTAT, N_target] deltaS;
}

transformed parameters {
  real rho_rec = tanh(z_rec);
  real rho_eq = tanh(z_eq);
  real rho_stat = tanh(z_stat);
}

model {
  // // prior distrributions for standard deviations
  phi_ss_sm ~ normal(0, 0.5);
  phi_ss_lm ~ normal(0, 0.5);
  phi_s2s ~ normal(0, 0.5);
  tau_sm ~ normal(0, 0.5);
  tau_lm ~ normal(0, 0.5);

  c0 ~ normal(0,0.5);

  z_eq ~ normal(pars_z_eq[1], pars_z_eq[2]);
  z_stat ~ normal(pars_z_stat[1], pars_z_stat[2]);
  z_rec ~ normal(pars_z_rec[1], pars_z_rec[2]);

  vector[N] phi_1 = phi_ss_sm[1] * M1_rec + phi_ss_lm[1] * M2_rec;
  vector[N] phi_2 = phi_ss_sm[2] * M1_rec + phi_ss_lm[2] * M2_rec;
  vector[NEQ] tau_1 = tau_sm[1] * M1_eq + tau_lm[1] * M2_eq;
  vector[NEQ] tau_2 = tau_sm[2] * M1_eq + tau_lm[2] * M2_eq;

  deltaB[:,1] ~ normal(0, tau_1);
  deltaS[:,1] ~ normal(0, phi_s2s[1]);

  // frst target
  vector[N] mu_1 = c0[1] + deltaB[eq, 1] + deltaS[stat,1];
  Y[:,1]  ~ normal(mu_1, phi_1);
  vector[N] deltaWS_1 = Y[:,1] - mu_1;

  // second target
  vector[NEQ] mu_deltaB = tau_2 ./ tau_1 * rho_eq .* deltaB[:,1];
  vector[NEQ] tau_cond = sqrt((1 - square(rho_eq)) * square(tau_2));
  deltaB[:,2] ~ normal(mu_deltaB, tau_cond);


  vector[NSTAT] mu_deltaS = phi_s2s[2] / phi_s2s[1] * rho_stat * deltaS[:,1];
  real phi_s2s_cond = sqrt((1 - square(rho_stat)) * square(phi_s2s[2]));
  deltaS[:,2] ~ normal(mu_deltaS, phi_s2s_cond);

  vector[N] mu_2 = c0[2] + deltaB[eq, 2] + deltaS[stat,2] +
                 phi_2 ./ phi_1 * rho_rec .* deltaWS_1;
  vector[N] phi_ss_cond = sqrt((1 - square(rho_rec)) * square(phi_2));

  Y[idx_obs_var2,2]  ~ normal(mu_2[idx_obs_var2], phi_ss_cond[idx_obs_var2]);
}

