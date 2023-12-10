#include<Rcpp.h>
using namespace Rcpp;
//' @title Based on EM algorithm and boot algorithm, the survival analysis of censored data was carried out
//' @description Based on EM algorithm and boot algorithm, the survival analysis of censored data was carried out
//' @param x the results of the first type of experiment
//' @param s cut-off values in the first type of experiment
//' @param y the results of the second type of experiment
//' @param u the results of the third type of experiment
//' @param v the results of the third type of experiment
//' @param lambda0 the initial value of the iteration
//' @param max.it the maximum number of iterations
//' @param lambda parameters of the survival model
//' @param B number of bootstrap samples
//' @param alpha the level of significance of the confidence interval
//' @return  a list with maximum likelihood estimates for parameters, full values for censored data, and estimated expectations, variances, and confidence intervals
//' @export
// [[Rcpp::export]]
List EMcpp(NumericVector x, double s, NumericVector y, NumericVector u, NumericVector v,
           double lambda0 = 0.1, int max_it = 1e4, double tol = 1e-6) {
  
  int m = x.length();
  int r = m - sum(x);
  int n = y.length();
  int k = u.length();
  
  NumericVector lambda(max_it);
  lambda[0] = lambda0;
  
  for (int i = 1; i < max_it; i++) {
    double n1 = m + n + k - 1;
    double n2 = s * exp(-s * lambda[i - 1]);
    double d2 = 1 - exp(-s * lambda[i - 1]);
    
    NumericVector n3 = (u * exp(-lambda[i - 1] * u) - v * exp(-lambda[i - 1] * v));
    NumericVector d3 = exp(-lambda[i - 1] * u) - exp(-lambda[i - 1] * v);
    
    double d1 = n * mean(y) + (m - r) * (s + 1 / lambda[i - 1]) + r * (1 / lambda[i - 1] - n2 / d2) + k / lambda[i - 1] + sum(n3 / d3);
    lambda[i] = n1 / d1;
  
  // Check convergence
  if (std::abs(lambda[i] - lambda[i - 1]) < tol) {
    return Rcpp::List::create(Rcpp::Named("lambda") = lambda[i], Rcpp::Named("iteration") = i + 1);
  }
  }
  
  // If convergence criteria are not met, return the result at the maximum iteration
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda[max_it - 1], Rcpp::Named("iteration") = max_it);
}


// [[Rcpp::export]]
NumericVector Supvaluecpp(NumericVector x, double s, NumericVector y, NumericVector u, NumericVector v, double lambda) {
  
  int m = x.length();
  int r = m - sum(x);
  int n = y.length();
  int k = u.length();
  
  NumericVector Y(m + n + k);
  
  NumericVector u1 = runif(r, 0, 1);
  NumericVector u2 = runif(m - r, 0, 1);
  
  for (int i = 0; i < r; i++) {
    Y[i] = -1 / lambda * log(1 - u1[i] * (1 - exp(-lambda * s)));
  }
  
  for (int i = 0; i < (m - r); i++) {
    Y[r + i] = s - log(u2[i]) / lambda;
  }
  
  for (int i = 0; i < n; i++) {
    Y[m + i] = y[i];
  }
  
  for (int i = 0; i < k; i++) {
    Y[m + n + i] = ((u[i] + 1 / lambda) * exp(-lambda * u[i]) - (v[i] + 1 / lambda) * exp(-lambda * v[i])) / (exp(-lambda * u[i]) - exp(-lambda * v[i]));
  }
  
  return Y;
}

List Bootstrapcpp(NumericVector x, int B, double alpha = 0.05) {
  
  // Calculate lambda hat
  double lambdahat = 1 / mean(x);
  
  // Vector to store bootstrap estimates
  NumericVector lambdastar(B);
  
  // Bootstrap resampling
  for(int b = 0; b < B; b++) {
    NumericVector xstar = sample(x, x.size(), true);
    lambdastar[b] = 1 / mean(xstar);
  }
  
  // Calculate bootstrap statistics
  double mean_boot = mean(lambdastar);
  double bias_boot = mean_boot - lambdahat;
  double se_boot = sd(lambdastar);
  
  // Confidence intervals
  //double normal_CI1 = lambdahat - Rcpp::qnorm(1 - 0.5 * alpha) * se_boot;
  //double normal_CI2 = lambdahat + Rcpp::qnorm(1 - 0.5 * alpha) * se_boot;
  
  //double basic_CI1 = 2 * lambdahat - Rcpp::quantile(lambdastar, 1 - 0.5 * alpha)[0];
  //double basic_CI2 = 2 * lambdahat - Rcpp::quantile(lambdastar, 0.5 * alpha)[0];
  
  //double percentile_CI1 = Rcpp::quantile(lambdastar, 0.5 * alpha)[0];
  //double percentile_CI2 = Rcpp::quantile(lambdastar, 1 - 0.5 * alpha)[0];
  
  // Create and return a List with named elements
  return List::create(Named("mean_boot", mean_boot),
                      Named("bias_boot", bias_boot),
                      Named("se_boot", se_boot));
}
