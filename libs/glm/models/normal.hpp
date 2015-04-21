#ifndef __NORMAL_H__
#define __NORMAL_H__

#include <glm/models/glm_model.hpp>
#include <glm/models/links/glm_link.hpp>

/**
 * Implements the normal linear regression model.
 */
class normal : public glm_model
{
public:
    /**
     * Constructor.
     *
     * @param link_name Name of the link function.
     */
    normal(const std::string &link_name);
    
    /**
     * Constructor.
     *
     * @param link A link function, this class will be responsible for freeing it.
     */
    normal(glm_link *link);

    /**
     * Destructor.
     */
    ~normal();

    /**
     * @see glm_model.valid_mu.
     */
    virtual bool valid_mu(const arma::vec &mu) const;

    /**
     * @see glm_model.compute_mu.
     */
    virtual arma::vec var(const arma::vec &mu) const;
    
    /**
     * @see glm_model.dispersion.
     */
    virtual double dispersion(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float k) const;

    /**
     * @see glm_model.likelihood.
     */
    virtual double likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float dispersion = 1.0) const;
};

#endif /* End of __NORMAL_H__ */
