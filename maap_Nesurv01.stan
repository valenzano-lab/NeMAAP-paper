data{
  int N; // Number of observations
  vector[N] x; // Predictor variable, Ne
  vector[N] y; // Outcome/response variable, median lifespan
}

transformed data{
  vector[N] logx; // Transformed data as the negative log of the population size
  logx = log(x); // Formula for the transformation logarithm
}

parameters{
  real alpha; // Parameter for slope 
  real beta; // Parameter for intercept
  real <lower=0> sigma; // Parameter for dispersion (standard deviation)
}

model{
  y ~ normal(alpha*logx + beta, sigma); // likelihood of the model
}

