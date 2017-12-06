
data {             
  int<lower=0> I;                          // # respondents
  int<lower=0> K;                          // # subpopulations
  vector[K] mu_beta;                       // fixed mean for beta priors
  vector<lower=0>[K] sigma_beta;           // fixed variance for beta priors
  int  y[I,K];                             // # known by person i in subpopulation k
}
parameters {
  vector[I] alpha;                         // log degree
  vector[K] beta;                          // log prevalence of group in population
  matrix<lower=0>[I,K] gam;                // relative propensity
  real mu_alpha;                           // prior mean for beta
  real<lower=0> sigma_alpha;               // prior scale for beta
}
model {
// priors
  alpha ~ normal(mu_alpha, sigma_alpha);   // identifies location and scale
  beta ~ normal(mu_beta, sigma_beta);      // hierarchical priors
  for (i in 1:I) for (k in 1:K) gam[i,k] ~ gamma(.25, .25);

// hyperpriors
  mu_alpha ~ normal(0, 10);               // weakly informative 
  sigma_alpha ~ inv_gamma(.001, .001);    // weakly informative 

  for (k in 1:K) {
    for (i in 1:I) {
      real xi_i_k;
      xi_i_k = gam[i,k] * exp(alpha[i] + beta[k])  ;
      y[i,k] ~ poisson(xi_i_k);
    }
  }
}


