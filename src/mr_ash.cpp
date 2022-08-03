#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void updatebetaj (const vec& xj, double wj, double& betaj, vec& r, vec& piold,
		  vec& pi, double sigma2, const vec& sa2, const vec& s2inv,
		  double& a1, double& a2, int j, int p, double epstol);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mr_ash_rcpp (const arma::mat& X, const arma::vec& y,
		  const arma::vec& w, const arma::vec& sa2,
		  arma::vec& pi, arma::vec& beta, arma::vec& r, 
		  double sigma2, const arma::uvec& o, int maxiter, 
		  int miniter, double convtol, double epstol, 
		  std::string method_q, bool updatepi, 
		  bool updatesigma, int verbose) {
  
  // ---------------------------------------------------------------------
  // DEFINE SIZES
  // ---------------------------------------------------------------------
  int n = X.n_rows;
  int p = X.n_cols;
  int K = sa2.n_elem;
  
  // ---------------------------------------------------------------------
  // PREDEFINE LOCAL VARIABLES
  // ---------------------------------------------------------------------
  vec varobj(maxiter);
  vec dbeta(maxiter);
  vec sigma2byiter(maxiter);
  vec w1(maxiter);
  unsigned int iter;
  int i                  = 0;
  
  double a1;
  double a2;
  vec piold;
  vec betaold;
  
  // ---------------------------------------------------------------------
  // PRECALCULATE
  // ---------------------------------------------------------------------
  mat S2inv        = 1 / outerAddition(1/sa2, w);
  S2inv.row(0).fill(epstol);
  
  // Repeat until convergence criterion is met, or until the maximum
  // number of iterations is reached.
  for (iter = 0; iter < maxiter; iter++) {
    
    // reset parameters
    a1                   = 0;
    a2                   = 0;
    piold                = pi;
    pi.fill(0);
    betaold              = beta;
    
    // ---------------------------------------------------------------------
    // RUN COORDINATE ASCENT UPDATES : INDEX 1 - INDEX P
    // ---------------------------------------------------------------------
    for (unsigned int j = 0; j < p; j++){
      
      updatebetaj(X.col(o(i)), w(o(i)), beta(o(i)), r, piold, pi, sigma2, sa2, S2inv.col(o(i)), a1, a2, o(i), p, epstol);
      i++;
    }
    
    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 1
    // ---------------------------------------------------------------------
    varobj(iter) = dot(r,r) - dot(square(beta), w) + a1;
    
    // ---------------------------------------------------------------------
    // UPDATE SIGMA2 IF REQUESTED
    // ---------------------------------------------------------------------
    if (updatesigma) {
      if (method_q == std::string("sigma_indep_q")) {
        sigma2 = varobj(iter) + p * (1.0 - pi(0)) * sigma2;
        sigma2 = sigma2 / (n + p * (1.0 - pi(0)));
      } else if (method_q == std::string("sigma_dep_q")) {
        sigma2 = varobj(iter) / n;
      }
    }
    
    if (updatepi)
      piold = pi;

    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 2
    // ---------------------------------------------------------------------
    varobj(iter) = varobj(iter) / sigma2 / 2.0 +
                   log(2.0 * M_PI * sigma2) / 2.0 * n -
                   dot(pi, log(piold + epstol)) * p + a2;
    
    for (unsigned int j = 1; j < K; j++)
      varobj(iter) += pi(j) * log(sa2(j)) * p / 2;
    
    if (!updatepi)
      pi = piold;
    
    // ---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // ---------------------------------------------------------------------
    dbeta(iter)        = norm(beta - betaold);
    sigma2byiter(iter) = sigma2;
    w1(iter)           = sum(pi > 1e-4);
    if (verbose == 1)
      Rprintf("+");
    else if (verbose == 2)
      Rprintf("%4d %+0.12e %0.2e %0.2e %3d\n",iter + 1,varobj(iter),
	      dbeta(iter),sqrt(sigma2),w1(iter));
    if (iter >= miniter - 1) {
      if (norm(betaold - beta) < convtol * p) {
        iter++;
        break;
      }
      if (iter > 0) {
        if (varobj(iter) > varobj(iter - 1))
          break;
      }
    }
  }
  
  if (verbose == 1)
    Rprintf("\n");
  if (verbose > 0)
    Rprintf("Mr.ASH terminated at iteration %d.\n", iter);

  // ---------------------------------------------------------------------
  // RETURN VALUES
  // ---------------------------------------------------------------------
  return List::create(Named("beta")   = beta,
                      Named("sigma2") = sigma2,
                      Named("pi")     = pi,
                      Named("iter")   = iter,
                      Named("varobj") = varobj,
		      Named("dbeta")  = dbeta,
		      Named("sigma2byiter") = sigma2byiter,
		      Named("w1")     = w1);
}

// TO DO: Explain here what this function does, and how to use it.
void updatebetaj (const vec& xj, double wj, double& betaj, vec& r, vec& piold,
		  vec& pi, double sigma2, const vec& sa2, const vec& s2inv,
		  double& a1, double& a2, int j, int p, double epstol) {
  
  // calculate b
  double bjwj = dot(r, xj) + betaj * wj;
  
  // update r first step
  r += xj * betaj; 
  
  // calculate muj
  vec muj = bjwj * s2inv;
  muj(0) = 0;
  
  // calculate phij
  vec phij = log(piold + epstol) - log(1 + sa2 * wj)/2 + muj * (bjwj / 2 / sigma2);
  phij = exp(phij - max(phij));
  phij = phij / sum(phij);
  
  // pinew
  pi += phij / p;
  
  // update betaj
  betaj = dot(phij, muj);
  
  // update r second step
  r += -xj * betaj;
  
  // precalculate for M-step
  a1 += bjwj * betaj;
  a2 += dot(phij, log(phij + epstol));
  phij(0) = 0;
  a2 += -dot(phij, log(s2inv)) / 2;
  
  return;
}

