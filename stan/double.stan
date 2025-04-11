data {
  int<lower=1> N;
  array[N] int D1; // phenomenon detected flag
  array[N] int D2; // phenomenon detected flag
  array[N] int R1; // num. records
  array[N] int R2; // num. records
  real<lower=0,upper=1> q_mu;
  real<lower=0> q_theta; 
}
parameters {
  real<lower=0,upper=1> q;
  real<lower=0,upper=1> rho1;
  real<lower=0,upper=1> rho2;
}
model {
  q ~ beta(q_mu * q_theta, (1-q_mu) * q_theta);
  rho1 ~ beta(1, 3);
  rho2 ~ beta(1, 3);
  for (i in 1:N) {
    if (R1[i] > 0 && R2[i] == 0) {
      if (D1[i] == 0) {
        target += log(q * (1 - rho1)^R1[i] + (1 - q));
      } else {
        target += log(q * (1 - (1 - rho1)^R1[i]));
      }
    } else if (R1[i] == 0 && R2[i] > 0) {
      if (D2[i] == 0) {
        target += log(q * (1 - rho2)^R2[i] + (1 - q));
      } else {
        target += log(q * (1 - (1 - rho2)^R2[i]));
      }
    } else if (R1[i] > 0 && R2[i] > 0) {
      if (D1[i] == 0 && D2[i] == 0) {
        target += log(q * (1 - rho1)^R1[i] * (1 - rho2)^R2[i] + (1 - q));
      } else if (D1[i] == 1 && D2[i] == 1) {
        target += log(q * (1 - (1 - rho1)^R1[i]) * (1 - (1 - rho2)^R2[i]));
      } else if (D1[i] == 1 && D2[i] == 0) {
        target += log(q * (1 - (1 - rho1)^R1[i]) * (1 - rho2)^R2[i]);
      } else if (D1[i] == 0 && D2[i] == 1) {
        target += log(q * (1 - rho1)^R1[i] * (1 - (1 - rho2)^R2[i]));
      }
    }
  }
}
