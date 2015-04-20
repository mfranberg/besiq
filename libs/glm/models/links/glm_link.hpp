#ifndef __GLM_LINK_H__
#define __GLM_LINK_H__

#include <armadillo>

class glm_link
{
public:
    /**
     * Constructor.
     *
     * @param name The name of link function.
     */
    glm_link(const std::string &name)
        : m_name( name )
    {
    }

    /**
     * Destructor.
     */
    virtual ~glm_link(){ };

    /**
     * Returns the name of this link function.
     *
     * @return the name of this link function.
     */
    std::string get_name() const
    {
        return m_name;
    }

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
     * The derivative of mu with respect to eta. It is often good to
     * use the relationship dmu/deta = (deta/dmu)^-1.
     *
     * @param mu The mean value parameter.
     *
     * @return The derivative of mu with respect to eta.
     */
    virtual arma::vec mu_eta(const arma::vec &mu) const = 0;

    /**
     * Compute the mean value parameter from the linearized parameter.
     *
     * @param eta The linearized parameter.
     *
     * @return The mean value parameter.
     */
    virtual arma::vec mu(const arma::vec &eta) const = 0;

private:
    std::string m_name;
};

/**
 * Creates a link function from the given name.
 *
 * The following links are avaiable:
 * - "identity" g(mu) = mu
 * - "log" g(mu) = log(mu)
 * - "logc" g(mu) = log(1-mu)
 * - "logit" g(mu) = log(mu/(1-mu))
 * - "odds" g(mu) = mu/(1-mu)
 *
 * @param link_name The name of the link function.
 *
 * @return The link if found, or null otherwise.
 */
glm_link *make_link(const std::string &link_name);

#endif /* End of __GLM_LINK_H__ */
