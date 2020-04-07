
data {
  
  int<lower=1> N;  // total number of observations 
  //int<lower=0> K;  // total number of predictor columns
                   // equivalent to length of random effects vector
  int<lower=1> G;  // number of hazard intervals
  //int<lower=1> J;  // number of patients
  //int<lower=1> nf; // number of fixed effects
  //int<lower=1> nr; // number of random effects
  //int<lower=1> indexf[nf]; // indexes from X of fixed effects
  //int<lower=1> indexr[nr]; // indexes from X of random effects
  int<lower=0> c[N];  // counts of event (model response)
  //matrix[N,K] x;   // design matrix, no intercept
  //vector[N] H;      // offset (not logged)
  matrix[N,G] xl;  // design matrix for interval effects
  //int<lower=1,upper=J> jj[N];  // patient indicator for each observation
  
}

transformed data {
  // Recommended QR decomposition of the data x for better posterior inference
  //matrix[N,K] Qast;
  //matrix[K,K] Rast;
  //matrix[K,K] Rast_inv;
  // Thin and scale the QR decomposition
  //Qast = qr_Q(x)[, 1:K] * sqrt(N - 1);
  //Rast = qr_R(x)[1:K, ] / sqrt(N - 1);
  //Rast_inv = inverse(Rast);
}

parameters {
  // random patient effects covariance matrix components
  //vector<lower=0>[nr] sigma;  // variances on random effects
  //corr_matrix[K] Omega;       // correlation matrix
  //cholesky_factor_corr[nr] Lcorr;

  vector[G] lambda;  // intercepts 
  //vector[nr] theta;  // pop. coefficients of random effects
  //vector[nf] thetaf;  // coefficients of fixed effects
  //vector[nr] theta_random[J];   // J x nr of patient level coefficients

}

transformed parameters {
  
}

model {
  
  // prior for the interval risks
  //lambda ~ gamma(0.001, 0.001);
  lambda ~ normal(0, 5);
  
  // prior for population effects
  //theta ~ normal(0, 5);
  //thetaf ~ normal(0, 5);
  
  //prior on sigma
  //sigma ~ normal(0, 5); // is prior sd too small?
  
  // Prior on L-corr matrix
  //Lcorr ~ lkj_corr_cholesky(2);
  
  // multivariate prior dist on random effects
  //beta_random ~ multi_normal(zero_beta, Tau);
  //theta_random ~ multi_normal_cholesky(theta, diag_pre_multiply(sigma, Lcorr));
  
  {  // vectorized format of sampling distribution
     // for increased efficiency
     vector[N] p;
     
     // link
     for(n in 1:N) {
       p[n] = inv_cloglog( xl[n]*lambda ); 
     }
     
     // sampling distribution centered at the prediction equation
     c ~ bernoulli(p);
    
  }
}

generated quantities {
  // Generate covariance matrix
  //matrix[nr,nr] Omega;
  //matrix[nr,nr] Sigma;
  //vector[N] log_lik;
  
  //Omega = multiply_lower_tri_self_transpose(Lcorr);
  //Sigma = quad_form_diag(Omega, sigma);
  
  // Pointwise log-like for WAIC
  //for(n in 1:N) {
  // log_lik[n] = bernoulli_lpmf( c[n] | inv_cloglog( xl[n]*lambda ) );
  //}
}
