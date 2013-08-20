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
         * Compute the weight vector that is used in one iteration
         * in the Iteratively reweighted least squares algorithm.
         *
         * @param X The design matrix.
         * @param y The observations.
         * @param b The beta coefficients.
         *
         * @return A weight vector.
         */
        virtual arma::vec compute_w(const arma::mat &X, const arma::vec &y, const arma::vec &b) const = 0;

        /**
         * Compute the adjusted dependent variates in the Iteratively reweighted
         * least squares algorithm.
         *
         * @param X The design matrix.
         * @param y The observations.
         * @param b The beta coefficients.
         *
         * @return The adjusted dependent variates.
         */
        virtual arma::vec compute_z(const arma::mat &X, const arma::vec &y, const arma::vec &b) const = 0;

        /**
         * Compute the log likelihood for the parameters.
         *
         * @param X The design matrix.
         * @param y The observations.
         * @param b The beta coefficients.
         *
         * @return The log likelihood.
         */
        virtual double likelihood(const arma::mat &X, const arma::vec &y, const arma::vec &b) const = 0;
};

#endif /* End of __GLM_MODEL_H__ */
