#ifndef __POWER_H__
#define __POWER_H__

#include <glm/models/links/glm_link.hpp>

/**
 * Implements the power link function (box-cox).
 */
class power_link : public glm_link
{
public:
    /**
     * Constructor.
     */
    power_link(float lambda);

    /**
     * @see glm_model.init_beta.
     */
    virtual arma::vec init_beta(const arma::mat &X, const arma::vec &y) const;

    /**
     * @see glm_model.mu.
     */
    virtual arma::vec mu(const arma::vec &eta) const;
    
    /**
     * @see glm_model.eta.
     */
    virtual arma::vec eta(const arma::vec &mu) const;

    /**
     * @see glm_model.mu_eta.
     */
    virtual arma::vec mu_eta(const arma::vec &mu) const;

private:
    /**
     * Choose how to transform.
     */
    float m_lambda;
};

#endif /* End of __POWER_H__ */
