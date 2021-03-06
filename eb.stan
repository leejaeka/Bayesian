
functions {
  vector rhs(real t, vector y, real k1, real k2, real gamma, real P, real u){
    vector[3] dydt;
    dydt[1] = (40*1e-4) - 2*y[1] - k1*y[2]*y[1];
    dydt[2] = P*k2*y[3] - u*y[2] - k1*y[2]*y[1];
    dydt[3] = k1*y[2]*y[1] - gamma*y[3] - k2*y[3];
    return dydt;
  }
}

data {
  int<lower=0> N;
  vector[N] y;
  real t[N];
  real<lower=0> D; // initial dose
  real<lower=0> k1;
  real<lower=0> k2;
  real<lower=0> gamma;
  real<lower=0> sigma;
  real<lower=0> P;
  real<lower=0> u;
}

transformed parameters{
  
  vector[N] mu;
  {
    vector[2] solution[N] = ode_rk45(rhs, [9600*1e-4, 0.1, 0.0]', 0.0, t, k1, k2, gamma, P, u);
    
    for(i in 1:N){
      mu[i] = solution[i,2];
    }
  }
}

model {
  y ~ normal(mu, sigma);
  sigma ~ normal(0, 1);
  c_0 ~ normal(0, 1);
}

