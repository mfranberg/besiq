#ifndef __LOG_H__
#define __LOG_H__

#include <glm/models/links/glm_link.hpp>

/**
 * Implements the multiplicative penetrance model.
 */
class log_link : public glm_link
{
public:
    /**
     * Constructor.
     */
    log_link();

    /**
     * @see glm_model.init_beta.
     */
    virtual arma::vec init_beta(const arma::mat &X, const arma::vec &y) const;

    /**
     * @see glm_model.mu.
     */
    virtual arma::vec mu(const arma::vec &eta) const;

    /**
     * @see glm_model.mu_eta.
     */
    virtual arma::vec mu_eta(const arma::vec &mu) const;

private:
};

#endif /* End of __LOG_H__ */
