#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ------ Stein features ------
// [[Rcpp::export]]
arma::mat efficient_getX(
    const arma::mat& samples,
    const arma::mat& derivatives,
    const arma::mat& combinations
) {
  if (combinations.min() < 0 ||
      arma::any(arma::vectorise(combinations != arma::floor(combinations))))
    Rcpp::stop("combinations must be non-negative integers");

  const uword N = samples.n_rows, d = samples.n_cols, K = combinations.n_cols;
  const uword D = (uword) combinations.max();

  arma::cube P(N, d, D + 1);
  P.slice(0).ones();
  for (uword p = 1; p <= D; ++p) P.slice(p) = P.slice(p - 1) % samples;

  arma::mat X(N, K, fill::zeros);
  std::vector<uword> idx, pw;
  idx.reserve(d); pw.reserve(d);
  arma::vec term(N), rest(N);

  for (uword k = 0; k < K; ++k) {
    idx.clear(); pw.clear();
    for (uword j = 0; j < d; ++j) {
      const uword a = (uword) combinations(j, k);
      if (a > 0) { idx.push_back(j); pw.push_back(a); }
    }
    const uword s = idx.size();
    if (s == 0) continue;

    term.zeros();
    for (uword r = 0; r < s; ++r) {
      rest.ones();
      for (uword q = 0; q < s; ++q)
        if (q != r) rest %= P.slice(pw[q]).col(idx[q]);

        const uword j = idx[r];
        const double a = pw[r];
        term += a * rest % P.slice(pw[r] - 1).col(j) % derivatives.col(j);
        if (pw[r] > 1)
          term += a * (a - 1.0) * rest % P.slice(pw[r] - 2).col(j);
    }
    X.col(k) = -0.5 * term;
  }

  return X;
}

// ------ OLS with LOO residuals ------
// [[Rcpp::export]]
Rcpp::List fast_ols_loo_cpp(const arma::mat& Y, const arma::mat& Z) {
  const uword n = Z.n_rows, p = Z.n_cols + 1, M = Y.n_cols;

  arma::mat X(n, p);
  X.col(0).ones();
  X.cols(1, p - 1) = Z;

  const double ridge = 1e-8;
  arma::mat Xaug = arma::join_cols(X, std::sqrt(ridge) * arma::eye(p, p));
  arma::mat Yaug = arma::join_cols(Y, arma::zeros(p, M));

  arma::mat Q, R;
  arma::qr_econ(Q, R, Xaug);
  arma::mat beta = arma::solve(arma::trimatu(R), Q.t() * Yaug);

  arma::mat E = Y - X * beta;

  arma::vec h = arma::sum(arma::square(Q.rows(0, n - 1)), 1);
  h = arma::clamp(h, 0.0, 1.0 - 1e-10);

  arma::mat loo_E = E.each_col() / (1.0 - h);

  return Rcpp::List::create(
    Rcpp::Named("intercept") = beta.row(0),
    Rcpp::Named("loo_residuals") = loo_E
  );
}
