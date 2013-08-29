#ifndef __GLM_MODEL_H__
#define __GLM_MODEL_H__

#include <armadillo>

/**
 * General base class for GLM models that will be solved with
 * the Iteratively reweighted least squares algorithm.
 */
class glm_model
{
    public:
        /**
         * Generate starting values for the beta coefficients.
         *
         * @param X The design matrix, including an intercept if desired.
         * @param y The outcome or observations.
         *
         * @return A vector containing the starting values for beta.
         */
        virtual arma::vec init_beta(const arma::mat &X, const arma::vec &y) const = 0;

        /**
         * Compute the mean value parameter from the linearized parameter.
         *
         * @param eta The linearized parameter.
         *
         * @return The mean value parameter.
         */
        virtual arma::vec mu(const arma::vec &eta) const = 0;
        
        /**
         * The derivative of mu with respect to eta. It is often good to
         * use the relationship dmu/deta = (deta/dmu)^-1.
         *
         * @param mu The mean value parameter.
         *
         * @return The derivative of mu with respect to eta.
         */
        virtual arma::vec mu_eta(const arma::vec &mu) const = 0;

        /**
         * The variance of each observation given the mean value
         * parameter.
         *
         * @param mu The mean value parameter.
         *
         * @return The variance of each observation.
         */
        virtual arma::vec var(const arma::vec &mu) const = 0;

        /**
         * Compute the log likelihood for the parameters.
         *
         * @param mu The mean value parameter.
         * @param y The observations.
         * @param missing Indiciates missing samples by 1 and not missing by 0.
         *
         * @return The log likelihood.
         */
        virtual double likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const = 0;
};

#endif /* End of __GLM_MODEL_H__ */
