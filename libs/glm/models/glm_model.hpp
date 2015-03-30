#ifndef __GLM_MODEL_H__
#define __GLM_MODEL_H__

#include <armadillo>

class glm_link;

/**
 * General base class for GLM models that will be solved with
 * the Iteratively reweighted least squares algorithm.
 */
class glm_model
{
    public:
        /**
         * Returns the link function used in this model.
         *
         * @return the link function used in this model.
         */
        virtual const glm_link &get_link() const = 0;

        /**
        * Checks that mu is in the allowed range. Should return false
        * if this is not the case, the irls algorithm will terminate
        * in this case.
        *
        * @param mu The mean value parameter.
        *
        * @return True if mu is in range, false otherwise.
        */
        virtual bool valid_mu(const arma::vec &mu) const = 0;
        
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
