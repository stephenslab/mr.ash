#include "misc.h"

using namespace arma;

// [[Rcpp::export]]
arma::uvec random_order (int p, int numiter) {
  uvec o(p * numiter);
  for (unsigned int i = 0 ; i < numiter; i++)
    o.subvec(i * p, (i+1) * p - 1) = randperm(p);
  return o;
}

mat outerAddition (const vec& a, const vec& b) {
  mat A(a.n_elem, b.n_elem);
  A.fill(0);
  A.each_row()          += b.t();
  A.each_col()          += a;
  return A;
}
