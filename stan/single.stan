data {
  int<lower=1> N;
  array[N] int D; // phenomenon detected flag
  array[N] int R; // num. records
  real<lower=0,upper=1> q_mu;
  real<lower=0> q_theta; 
}
parameters {
  real<lower=0,upper=1> q;
  real<lower=0,upper=1> rho;
}
model {
  q ~ beta(q_mu * q_theta, (1-q_mu) * q_theta);
  rho ~ beta(1, 3);
  for (i in 1:N) {
    if (R[i] > 0) {
      if (D[i] == 0) {
        target += log(q * (1 - rho)^R[i] + (1 - q));
      } else {
        target += log(q * (1 - (1 - rho)^R[i]));
      }
    }
  }
}
