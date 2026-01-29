#include <TMB.hpp>

// numerically stable log(1 - Phi(z))
template<class Type>
Type log1mPhi(Type z) {
  // For very large z, use tail approximation
  if(z > Type(6.0)) {
    return dnorm(z, Type(0.0), Type(1.0), true) - log(z);
  } else {
    return log(Type(1.0) - pnorm(z));
  }
}

// truncated normal log-density
template<class Type>
Type dtruncnorm(Type y, Type mu, Type sigma, Type L, bool give_log = false) {

  Type z  = (y - mu) / sigma;
  Type zL = (L - mu) / sigma;

  // log φ(z)
  Type logphi = dnorm(z, Type(0.0), Type(1.0), true);

  // log survival
  Type logS = log1mPhi(zL);

  Type out = logphi - logS - log(sigma);

  return give_log ? out : exp(out);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_VECTOR(y);         // outcome
  DATA_MATRIX(X);         // fixed effects design
  DATA_MATRIX(Z);         // random effects design
  DATA_IVECTOR(id);       // subject index (0-based)
  DATA_VECTOR(L);         // left truncation times

  int n = y.size();
  int p = X.cols();
  int q = Z.cols();
  int n_subjects = id.maxCoeff() + 1;

  // PARAMETERS
  PARAMETER_VECTOR(beta);          // fixed effects
  PARAMETER(log_sigma);            // residual SD
  PARAMETER_VECTOR(b_raw);         // stacked random effects (standard normal)
  PARAMETER_VECTOR(log_sd_b);      // vector of SDs for random effects

  Type sigma = exp(log_sigma);

  // Random effects structure:
  // b_raw is length n_subjects * q
  // b = diag(sd_b) * b_raw_subject
  vector<Type> sd_b = exp(log_sd_b);

  // Negative log-likelihood
  Type nll = 0.0;

  // Random effects contribution
  for(int i = 0; i < n_subjects * q; i++) {
    nll -= dnorm(b_raw(i), Type(0.0), Type(1.0), true);
  }

  // Observation likelihood
  for(int i = 0; i < n; i++) {

    int subj = id(i);
    // extract subject-level random effects
    vector<Type> b_i(q);
    for(int k = 0; k < q; k++) {
      b_i[k] = b_raw[subj*q + k] * sd_b[k];
    }

    // linear predictor
    Type mu = Type(0.0);
    for(int j = 0; j < p; j++) mu += X(i,j) * beta(j);
    for(int j = 0; j < q; j++) mu += Z(i,j) * b_i(j);

    // left truncation likelihood contribution
    nll -= dtruncnorm(y(i), mu, sigma, L(i), true);
  }

  return nll;
}
