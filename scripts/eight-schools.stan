data{

  int<lower = 0> J; //the number of schools

  real y[J]; // estimated treatment effect in each school

  real<lower=0> sigma[J]; // standard errors of effect estimates

  int test_position[3];//
  //
  vector[3] thing;


}

parameters{

real mu; // the average treatment effect

real<lower = 0> tau; // scalar for individual school effects

vector[J] eta; // vector of individual school effects

}

transformed parameters{

  vector[J] theta;

  vector[J*2] test;

  theta = mu + tau * eta;

   test[1:J] = theta;
  //
  test[(J + 1):(2*J)] = theta * 2;

  test[test_position] = thing;

  print(test[test_position])
  // test[3] = 1;

  // test[1] = 1;

  //test = append_row(test, theta);

  // print(test)

}

model{

target += normal_lpdf(eta | 0,1);

target += normal_lpdf(y | theta, sigma);

}
