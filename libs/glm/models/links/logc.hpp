#ifndef __LOGC_H__
#define __LOGC_H__

#include <glm/models/links/glm_link.hpp>

/**
 * Implements the logcomplement log(1 - p) model.
 */
class logc_link : public glm_link
{
public:
    /**
     * Constructor.
     */
    logc_link();

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

#endif /* End of __LOGC_H__ */
