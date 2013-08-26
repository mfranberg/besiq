#ifndef __DIRICHLET_H__
#define __DIRICHLET_H__

#include <armadillo>

#include <cmath>
#include <numeric>
#include <vector>

double dirmult(const arma::vec &x, const arma::vec &alpha);
double ldirmult(const arma::vec &x, const arma::vec &alpha);

#endif /* End of __DIRICHLET_H__ */
