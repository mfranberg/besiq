#ifndef __GLM_MODEL_H__
#define __GLM_MODEL_H__

#include <armadillo>

#include <glm/models/links/glm_link.hpp>

/**
 * General base class for GLM models that will be solved with
 * the Iteratively reweighted least squares algorithm.
 */
class glm_model
{
public:
    /**
     * Constructor.
     *
     * @param model_name Name of the model.
     * @param link_name Name of the link function.
     */
    glm_model(const std::string &model_name, const std::string &link_name)
        : m_model( model_name )
    {
        m_link = make_link( link_name );
    }

    /**
     * Constructor.
     *
     * @param model_name Name of the model.
     * @param link A link function, this class will be responsible for freeing it.
     */
    glm_model(const std::string &model_name, glm_link *link)
        : m_model( model_name )
    {
        m_link = link;
    }

    /**
     * Destructor.
     */
    virtual ~glm_model()
    {
        delete m_link;
    };

    /**
     * Returns the name of the model.
     *
     * @return the name of the model.
     */
    std::string get_name() const
    {
        return m_model;
    }

    /**
     * Returns the link function used in this model.
     *
     * @return the link function used in this model.
     */
    const glm_link &get_link() const
    {
        return *m_link;
    }

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
     * Estimate the dispersion of the model.
     *
     * @param mu The mean value parameter.
     * @param y The observations.
     * @param missing Indicates missing samples by 1 and not missing by 0.
     * @param k The number of estimated betas.
     *
     * @return The estimated dispersion.
     */
    virtual double dispersion(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float k) const = 0;

    /**
     * Compute the log likelihood for the parameters.
     *
     * @param mu The mean value parameter.
     * @param y The observations.
     * @param missing Indicates missing samples by 1 and not missing by 0.
     * @param dispersion Estimated dispersion (only used for final likelihood).
     *
     * @return The log likelihood.
     */
    virtual double likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float dispersion = 1.0) const = 0;
    
private:
    /**
     * Name of the model.
     */
    std::string m_model;

    /**
     * The link function.
     */
    glm_link *m_link;        
};

#endif /* End of __GLM_MODEL_H__ */
