#ifndef __PENETRANCE_MULTIPLICATIVE_H__
#define __PENETRANCE_MULTIPLICATIVE_H__

#include <glm/models/glm_model.hpp>

/**
 * Implements the multiplicative penetrance model.
 */
class penetrance_multiplicative : public glm_model
{
public:
    /**
     * @see glm_model.init_beta.
     */
    virtual arma::vec init_beta(const arma::mat &X, const arma::vec &y) const;

    /**
     * @see glm_model.mu.
     */
    virtual arma::vec mu(const arma::vec &eta) const;

    /**
     * @see glm_model.valid_mu.
     */
    virtual bool valid_mu(const arma::vec &mu) const;

    /**
     * @see glm_model.mu_eta.
     */
    virtual arma::vec mu_eta(const arma::vec &mu) const;

    /**
     * @see glm_model.compute_mu.
     */
    virtual arma::vec var(const arma::vec &mu) const;

    /**
     * @see glm_model.likelihood.
     */
    virtual double likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const;

private:
};

#endif /* End of __PENETRANCE_MULTIPLICATIVE_H__ */
