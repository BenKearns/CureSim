// The following has the code for all four dynamic models used.
// For implementation, each model has its own '.stan' file which is compiled prior to running.

// ******************** Cure fraction models ********************

// Damped trend cure model

data {
    int<lower=1> T;  // Time points
    vector[T-1] tau;   // Width between time-points (used for betas)
    vector[T] width;  // Widths used for gen pop ests
    int y[T];        // Events
    vector[T] n;     // At risk
    vector[T] Haz_pop;  // Gen pop hazard
    vector[T] Srv_pop;  // Gen pop survival
    real Max_cure;      // Max cure %
}

parameters {
    real beta_01;             // Initial coeff1
    real beta_02;             // Initial coeff2
    real<lower=0> Z;          // Variance coeff2
    vector[T-1] zeta_tilde;     // Tranformation of zeta2 (as in 8-schools example)
    real<lower=0, upper = Max_cure> cure_p;   // Cure proportion
    real<lower=0.8, upper=0.999> phi;
}

transformed parameters {
    vector[T] beta1;         // State 1
    vector[T-1] beta2;         // State 2
    vector[T] Srv_dis;       // Disease-specific survival
    vector[T] Srv_mod;      // Model estimate of overall survival
    vector[T] Haz_mod;      // Model estimate of overall hazard
    vector[T] cum_y;         // Used to calc Srv_dis
    
    { // Don't want to save this
      vector[T-1] zeta2;         // Innovations

      zeta2 = sqrt(Z) * zeta_tilde;
      beta1[1] = beta_01;
      beta2[1] = beta_02 + zeta2[1];
      for (t in 2:T-1) {
        beta1[t] = beta1[t-1] + beta2[t-1] * tau[t-1];
        beta2[t] = beta2[t-1] + zeta2[t];
      }
      beta1[T] = beta1[T-1] + beta2[T-1] * tau[T-1];
      cum_y = cumulative_sum(exp(beta1) .* width);
      Srv_dis = exp(-cum_y);
      Srv_mod = Srv_pop * cure_p + Srv_dis * (1 - cure_p);
      Haz_mod = (Haz_pop .* Srv_pop * cure_p + exp(beta1) .* Srv_dis * (1 - cure_p)) ./
                              (Srv_pop * cure_p + Srv_dis * (1 - cure_p));
    }
}

model {
  Z ~ inv_gamma(1, 0.005);
  zeta_tilde ~ normal(0, 1);
  y ~ poisson(Haz_mod .* n);
}

generated quantities{
  real level;
  real trend;
  real y_cum;

  level = beta1[T];
  trend = beta2[T-1];
  y_cum = cum_y[T];
}


// Local trend cure model

data {
    int<lower=1> T;  // Time points
    vector[T-1] tau;   // Width between time-points (used for betas)
    vector[T] width;  // Widths used for gen pop ests
    int y[T];        // Events
    vector[T] n;     // At risk
    vector[T] Haz_pop;  // Gen pop hazard
    vector[T] Srv_pop;  // Gen pop survival
    real Max_cure;      // Max cure %
}

parameters {
    real beta_01;             // Initial coeff1
    real beta_02;             // Initial coeff2
    real<lower=0> Z;          // Variance coeff2
    vector[T-1] zeta_tilde;     // Tranformation of zeta2 (as in 8-schools example)
    real<lower=0, upper = Max_cure> cure_p;   // Cure proportion
}

transformed parameters {
    vector[T] beta1;         // State 1
    vector[T-1] beta2;         // State 2
    vector[T] Srv_dis;       // Disease-specific survival
    vector[T] Srv_mod;      // Model estimate of overall survival
    vector[T] Haz_mod;      // Model estimate of overall hazard
    vector[T] cum_y;         // Used to calc Srv_dis
    
    { // Don't want to save this
      vector[T-1] zeta2;         // Innovations

      zeta2 = sqrt(Z) * zeta_tilde;
      beta1[1] = beta_01;
      beta2[1] = beta_02 + zeta2[1];
      for (t in 2:T-1) {
        beta1[t] = beta1[t-1] + beta2[t-1] * tau[t-1];
        beta2[t] = beta2[t-1] + zeta2[t];
      }
      beta1[T] = beta1[T-1] + beta2[T-1] * tau[T-1];
      cum_y = cumulative_sum(exp(beta1) .* width);
      Srv_dis = exp(-cum_y);
      Srv_mod = Srv_pop * cure_p + Srv_dis * (1 - cure_p);
      Haz_mod = (Haz_pop .* Srv_pop * cure_p + exp(beta1) .* Srv_dis * (1 - cure_p)) ./
                              (Srv_pop * cure_p + Srv_dis * (1 - cure_p));
    }
}

model {
  Z ~ inv_gamma(1, 0.005);
  zeta_tilde ~ normal(0, 1);
  y ~ poisson(Haz_mod .* n);
}

generated quantities{
  real level;
  real trend;
  real y_cum;

  level = beta1[T];
  trend = beta2[T-1];
  y_cum = cum_y[T];
}


// ******************** Non-cure models ********************

// Damped trend model

data {
    int<lower=1> T;  // Time points
    vector[T-1] tau;   // Width between time-points
    int y[T];        // Events
    vector[T] n;     // At risk
}

parameters {
    real beta_01;              // Initial coeff1
    real<lower=0> Z;          // Variance coeff2
    real beta_02;              // Initial coeff2
    vector[T-1] zeta_tilde;      // Tranformation of zeta2 (as in 8-schools example)
    real<lower=0.7, upper=0.999> phi;  // Whether or not we know phi values (don't let = 1 as messes up extrap calcs)
}

transformed parameters {
    vector[T] beta1;         // State 1
    vector[T-1] beta2;         // State 2
    { // Don't want to save this
      vector[T-1] zeta2;         // Innovations

      zeta2 = sqrt(Z) * zeta_tilde;
      beta1[1] = beta_01;
      beta2[1] = beta_02 + zeta2[1];
      for (t in 2:T-1) {
        beta1[t] = beta1[t-1] + beta2[t-1] * phi * tau[t-1];
        beta2[t] = beta2[t-1] * phi + zeta2[t];
      }
      beta1[T] = beta1[T-1] + beta2[T-1] * phi * tau[T-1];
    }
}

model {
  Z ~ inv_gamma(1, 0.005);
  zeta_tilde ~ normal(0, 1);
  y ~ poisson(exp(beta1) .* n);
}

generated quantities{
  real level;
  real trend;

  level = beta1[T];
  trend = beta2[T-1];
}


// Local trend model

data {
    int<lower=1> T;  // Time points
    vector[T-1] tau;   // Width between time-points
    int y[T];        // Events
    vector[T] n;     // At risk
}

parameters {
    real beta_01;              // Initial coeff1
    real<lower=0> Z;          // Variance coeff2
    real beta_02;              // Initial coeff2
    vector[T-1] zeta_tilde;      // Tranformation of zeta2 (as in 8-schools example)
}

transformed parameters {
    vector[T] beta1;         // State 1
    vector[T-1] beta2;         // State 2
    { // Don't want to save this
      vector[T-1] zeta2;         // Innovations

      zeta2 = sqrt(Z) * zeta_tilde;
      beta1[1] = beta_01;
      beta2[1] = beta_02 + zeta2[1];
      for (t in 2:T-1) {
        beta1[t] = beta1[t-1] + beta2[t-1] * tau[t-1];
        beta2[t] = beta2[t-1] + zeta2[t];
      }
      beta1[T] = beta1[T-1] + beta2[T-1] * tau[T-1];
    }
}

model {
  Z ~ inv_gamma(1, 0.005);
  zeta_tilde ~ normal(0, 1);
  y ~ poisson(exp(beta1) .* n);
}

generated quantities{
  real level;
  real trend;

  level = beta1[T];
  trend = beta2[T-1];
}

