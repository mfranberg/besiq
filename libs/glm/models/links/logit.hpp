#ifndef __LOGIT_H__
#define __LOGIT_H__

#include <glm/models/links/glm_link.hpp>

/**
 * Implements the binomial or logistic regression model.
 */
class logit_link : public glm_link
{
public:
    /**
     * Constructor.
     */
    logit_link();

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
};

#endif /* End of __LOGIT_H__ */
