#ifndef __BETA_H__
#define __BETA_H__

#include <armadillo>

#include <vector>

/**
 * Estimates parameters for the beta distribution using the
 * method of moments. The method of moments estimates will have
 * a higher variance than maximum likelihood, but is faster and
 * simpler.
 *
 * @param samples A list of samples to fit a beta distribution to.
 *
 * @return The estimated parameters, if they could not be estimated,
 *         the parameters 1,1 will be returned. 
 */
arma::vec mom_beta(const arma::vec &samples);

#endif /* End of __BETA_H__ */
