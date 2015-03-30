#ifndef __BINOMIAL_H__
#define __BINOMIAL_H__

#include <glm/models/glm_model.hpp>
#include <glm/models/links/glm_link.hpp>

/**
 * Implements the binomial or logistic regression model.
 */
class binomial : public glm_model
{
public:
    /**
     * Constructor.
     *
     * @param link_name Name of the link function.
     */
    binomial(const std::string &link_name);

    /**
     * Destructor.
     */
    ~binomial();

    /**
     * @see glm_model.get_link.
     */
    virtual const glm_link &get_link() const;

    /**
     * @see glm_model.valid_mu.
     */
    virtual bool valid_mu(const arma::vec &mu) const;

    /**
     * @see glm_model.compute_mu.
     */
    virtual arma::vec var(const arma::vec &mu) const;

    /**
     * @see glm_model.likelihood.
     */
    virtual double likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const;

private:
    /**
     * The link function.
     */
    glm_link *m_link;
};

#endif /* End of __BINOMIAL_H__ */
