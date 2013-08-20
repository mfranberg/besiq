#ifndef __BINOMIAL_H__
#define __BINOMIAL_H__

#include <models/glm_model.hpp>

/**
 * Implements the binomial or logistic regression model.
 */
class binomial : public glm_model
{
public:
    /**
     * @see glm_model.init_beta.
     */
    virtual arma::vec init_beta(const arma::mat &X, const arma::vec &y) const;

    /**
     * @see glm_model.init_beta.
     */
    virtual arma::vec compute_w(const arma::mat &X, const arma::vec &y, const arma::vec &b) const;

    /**
     * @see glm_model.init_beta.
     */
    virtual arma::vec compute_z(const arma::mat &X, const arma::vec &y, const arma::vec &b) const;

    /**
     * @see glm_model.init_beta.
     */
    virtual double likelihood(const arma::mat &X, const arma::vec &y, const arma::vec &b) const;

private:
};

#endif /* End of __BINOMIAL_H__ */
