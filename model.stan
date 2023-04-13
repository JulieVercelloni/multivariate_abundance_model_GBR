functions {
/** inverse ilr transformation of a vector x, using the inverse of the transpose of the V matrix of the ilr (tVinv)
*/
  vector ilrinv(matrix tVinv, vector x, int ntaxa) {
         vector[ntaxa] z;
         vector[ntaxa] y;
         z = exp(tVinv * x);
         y = z / sum(z);
         return y;
	}
}

data {
  int<lower = 0> ntaxa; // number of taxa
  int<lower = 0> nstills; // number of stills 
  matrix[ntaxa, ntaxa - 1] tVinv; //back-transformation matrix for ilr transformation
  int counts[nstills, ntaxa]; //observed counts
  vector[nstills] cyclone; //cyclone in binary (centered)
  vector[nstills] bleach; //bleaching in binary (centered)
  vector[nstills] both; //interaction cyclone and bleaching in binary (centered)
}

transformed data {
  int<lower = 1> s;
  s = ntaxa - 1;
}

parameters {
  vector[s] beta0; //intercept 
  vector[s] beta1; //cyclone effect
  vector[s] beta2; //bleaching effect
  vector[s] beta3; //cyclone and bleaching effect
  vector[s] z[nstills]; //transform into predicted logratio coordinates
  cholesky_factor_corr[s] LOmega; //Cholesky factor of prior correlation
  vector<lower=0>[s] tau; //prior scale on covariances
}

transformed parameters {
  cholesky_factor_cov[s] LSigma;
  vector[s] x[nstills]; //predicted logratio coordinates
  vector[ntaxa] rho[nstills]; //predicted relative abundances

  LSigma = diag_pre_multiply(tau, LOmega);
  for(i in 1:nstills){
    x[i] = beta0 + beta1 * cyclone[i] +beta2 * bleach[i] +beta3 * both[i] + LSigma * z[i];
    rho[i] = ilrinv(tVinv, x[i], ntaxa);
  } 

}

model {
  for(i in 1:nstills) {
    counts[i] ~ multinomial(rho[i]); // observation model
    z[i] ~ normal(0, 1);
  }
  tau ~ cauchy(0, 2.5);
  LOmega ~ lkj_corr_cholesky(2);
  beta0 ~ cauchy(0, 2.5);
  beta1 ~ cauchy(0, 2.5);
  beta2 ~ cauchy(0, 2.5);
  beta3 ~ cauchy(0, 2.5);
}

generated quantities{
  vector [nstills] log_lik;
  for (j in 1:nstills)
    log_lik[j] = multinomial_lpmf(counts[j] | rho[j]); //log likelihood for WAIC
}
