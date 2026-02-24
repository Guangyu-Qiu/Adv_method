data {
  int<lower=1> N; // number of rows (observation)
  vector[N] y; // outcome

  vector[N] drugs_0;
  vector[N] age_0;
  vector[N] race_0;
  vector[N] edu_0;
  vector[N] bmi_overweight;
  vector[N] bmi_underweight;
  vector[N] smoke_0;
}

parameters {
  real beta_0;
  real beta_drugs_0;
  real beta_age_0;
  real beta_race_0;
  real beta_edu_0;
  real beta_bmi_overweight;
  real beta_bmi_underweight;
  real beta_smoke_0;

  real<lower=0> sigma; // half-normal
}

model {

  // Non-informative priors
  beta_0 ~ normal(0, 100);

  beta_drugs_0 ~ normal(0, 100);
  beta_age_0 ~ normal(0, 100);
  beta_race_0 ~ normal(0, 100);
  beta_edu_0 ~ normal(0, 100);
  beta_bmi_overweight ~ normal(0, 100);
  beta_bmi_underweight ~ normal(0, 100);
  beta_smoke_0 ~ normal(0, 100);

  sigma ~ normal(0, 100); 

  for (n in 1:N) {
    y[n] ~ normal(
      beta_0
      + beta_drugs_0 * drugs_0[n]
      + beta_age_0 * age_0[n]
      + beta_race_0 * race_0[n]
      + beta_edu_0 * edu_0[n]
      + beta_bmi_overweight * bmi_overweight[n]
      + beta_bmi_underweight * bmi_underweight[n]
      + beta_smoke_0 * smoke_0[n],
      sigma
    );
  }
}