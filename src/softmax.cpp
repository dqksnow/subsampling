#include <RcppArmadillo.h>
// #include <cmath>
// #include <chrono>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



// [[Rcpp::export]]
arma::mat softmax_dL_sq_cpp(const arma::mat& X, const arma::mat& Y_matrix,
                            const arma::mat& P, const arma::vec& p, int K, int d,
                            double scale) {
  arma::mat S = Y_matrix - P;
  arma::mat dL(K * d, K * d);
  arma::mat X_t = X.t();
  arma::vec p_sq = arma::square(p);
  for (int i = 0; i < K; ++i) {
    for (int j = i; j < K; ++j) {
      arma::mat XSij = X.each_col() % (S.col(i) % S.col(j) / p_sq);
      arma::mat block = X_t * XSij;
      dL.submat(i * d, j * d, (i + 1) * d - 1, (j + 1) * d - 1) = block;
      if (i != j) {
        dL.submat(j * d, i * d, (j + 1) * d - 1, (i + 1) * d - 1) = block;
      }
    }
  }

  dL /= scale;

  return dL;
}

// [[Rcpp::export]]
arma::mat softmax_ddL_cpp(const arma::mat& X, const arma::mat& P,
                          const arma::vec& p, int K, int d, double scale) {
  arma::mat ddL(K * d, K * d);
  arma::mat X_t = X.t();
  for (int i = 0; i < K; ++i) {
    for (int j = i; j < K; ++j) {
      arma::mat XPj(d, d);
      if (i == j) {
        XPj = (X.each_col() % ((P.col(i) - P.col(i) % P.col(j)) / p));
      } else {
        XPj = (X.each_col() % ((-P.col(i) % P.col(j)) / p));
      }
      arma::mat block = X_t * XPj;
      ddL.submat(i * d, j * d, (i + 1) * d - 1, (j + 1) * d - 1) = block;
      if (i != j) {
        ddL.submat(j * d, i * d, (j + 1) * d - 1, (i + 1) * d - 1) = block;
      }
    }
  }
  ddL /= scale;
  return ddL;
}

// [[Rcpp::export]]
arma::mat softmax_Omega_cpp(const arma::mat& X, const arma::mat& P1,
                            const arma::vec& p, int K, int d, double scale) {
  arma::mat Omega = arma::zeros<arma::mat>(K * d, K * d);
  arma::vec P_sq = arma::sum(arma::square(P1), 1); // N * 1
  arma::mat P0 = P1.cols(1, P1.n_cols - 1); // exclude 1st column
  arma::mat X_t = X.t();

  for (int i = 0; i < K; ++i) {
    for (int j = i; j < K; ++j) {
      arma::mat XP;
      if (j == i) {
        XP = X.each_col() % ((P_sq + 1 - 2 * P0.col(i)) %
          arma::square(P0.col(i)) / p);
      } else {
        XP = X.each_col() % (((P_sq - P0.col(i) - P0.col(j)) %
          (P0.col(i) % P0.col(j))) / p);
      }
      arma::mat block = X_t * XP;
      Omega.submat(i * d, j * d, (i + 1) * d - 1, (j + 1) * d - 1) = block;
      if (i != j) {
        Omega.submat(j * d, i * d, (j + 1) * d - 1, (i + 1) * d - 1) = block;
      }
    }
  }
  Omega /= scale;

  return Omega;
}

