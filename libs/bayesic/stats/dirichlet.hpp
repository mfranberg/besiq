#ifndef __DIRICHLET_H__
#define __DIRICHLET_H__

#include <armadillo>

#include <cmath>
#include <numeric>
#include <vector>

/**
 * Computes the dirichlet multinomial probability of a vector
 * x with prior parameter alpha.
 * 
 * @param x The observations.
 * @param alpha The prior parameters of the dirichlet density.
 *
 * @return The posterior probability of x.
 */
double dirmult(const arma::vec &x, const arma::vec &alpha);

/**
 * Computes the dirichlet multinomial log probability of a vector
 * x with prior parameter alpha.
 * 
 * @param x The observations.
 * @param alpha The prior parameters of the dirichlet density.
 *
 * @return The log posterior probability of x.
 */
double ldirmult(const arma::vec &x, const arma::vec &alpha);

#endif /* End of __DIRICHLET_H__ */
