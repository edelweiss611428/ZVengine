data {
  int<lower=1> D;
  array[D] int<lower=0> data_counts;
  vector[D] prior_mean;
  matrix[D, D] L_prior_cov;
  real<lower=0> area;
}

parameters {
  vector[D] x;
}

model {
  x ~ multi_normal_cholesky(
    prior_mean,
    L_prior_cov
  );

  data_counts ~ poisson(
    area * exp(x)
  );
}
