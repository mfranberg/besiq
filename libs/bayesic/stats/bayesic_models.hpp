#ifndef __BAYESIC_MODELS_H__
#define __BAYESIC_MODELS_H__

#include <armadillo>

#include <plink_file.hpp>
#include <stats/dirichlet.hpp>
#include <stats/log_scale.hpp>
#include <stats/snp_count.hpp>

class model
{
public:
    /**
     * Returns the prior of the model.
     */
    virtual log_double prior( ) = 0;

    /**
     * Returns the likelihood of the data.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight) = 0;
};

class saturated
: public model
{
public:
    saturated(log_double prior);
    virtual log_double prior( );

    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

private:
    log_double m_prior;
};

class null
: public model
{
public:
    null(log_double prior);
    virtual log_double prior( );

    static arma::vec compute_alpha(const snp_row &row1, const snp_row &row2);

    static log_double snp_prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
    static log_double pheno_prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

private:
    log_double m_prior;
};

class ld_assoc
: public model
{
public:
    ld_assoc(log_double prior, bool is_first);
    virtual log_double prior( );

    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

private:
    log_double m_prior;
    bool m_is_first;
};

#endif /* End of __BAYESIC_MODELS_H__ */
