#ifndef __BAYESIC_MODELS_H__
#define __BAYESIC_MODELS_H__

#include <armadillo>

#include <plink/snp_row.hpp>
#include <bayesic/stats/dirichlet.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/stats/snp_count.hpp>

/**
 * This class represents a general model that can compute a
 * model likelihood for two snps.
 */
class model
{
public:
    /**
     * Constructor.
     *
     * @param prior The prior for the model.
     * @param alpha The prior parameters.
     */
    model(log_double prior, const arma::vec &alpha)
    : m_prior( prior ),
      m_alpha( alpha )
    {
    }

    /**
     * Destructor.
     */
    virtual ~model(){ };
    
    /**
     * Returns the prior of the model.
     *
     * @return The prior probability.
     */
    log_double prior( )
    {
        return m_prior;
    }

    /**
     * Returns the prior parameters.
     *
     * @return the prior parameters.
     */
    arma::vec get_alpha()
    {
        return m_alpha;
    }

    /**
     * Computes the likelihood of the phenotype under this
     * model.
     *
     * @param row1 The first snp.
     * @param row2 The §second snp.
     * @param phenotype The phenotype, discrete 0.0 and 1.0.
     * @param weight A weight for each sample, this will be used instead of 1.0 as a count.
     *
     * @return The likelihood of the snps.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight) = 0;
    
private:
    /**
     * The prior probability for the model.
     */
    log_double m_prior;

    /**
     * The beta prior parameters.
     */
    arma::vec m_alpha;
};

/**
 * This model represents a saturated model in the sense
 * that each cell in the penetrance table has it is own
 * parameter.
 */
class saturated
: public model
{
public:
    /**
     * Constructor.
     *
     * @param prior The prior probability for the model.
     */
    saturated(log_double prior, const arma::vec &alpha);
    
    /**
     * @see model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
};

/**
 * This class represents a model where no snp is associated with
 * the phenotype. The snps may however be in ld.
 */
class null
: public model
{
public:
    /**
     * Constructor.
     *
     * @param prior The prior probability for the model.
     */
    null(log_double prior, const arma::vec &alpha);

    /**
     * Computes a sensible prior alpha based on the minor allele frequencies
     * of both SNPs.
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     *
     * @return A vector of alphas that reflects the estiamted mafs.
     */
    static arma::vec compute_alpha(const snp_row &row1, const snp_row &row2);

    /**
     * Computes the probability of the snps under the null model.
     *
     * Note: Convenient to use in other models, as they always compute the
     *       null or prior probability of the snps.
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     * @param phenotype The phenotype, discrete 0.0 and 1.0.
     * @param weight A weight for each sample, this will be used instead of 1.0 as a count.
     *
     * @return the probability of the snps under the null model.
     */
    static log_double snp_prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
    
    /**
     * Computes the probability of the phenotype under the null model.
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     * @param phenotype The phenotype, discrete 0.0 and 1.0.
     * @param weight A weight for each sample, this will be used instead of 1.0 as a count.
     *
     * @return the probability of the phenotype under the null model.
     */
    static log_double pheno_prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

    /**
     * @see model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
};

/**
 * This class represents a model where one snp is associated with
 * the phenotype, and the other possibly in ld with the first.
 */
class ld_assoc
: public model
{
public:
    /**
     * Constructor.
     *
     * @param prior The prior probability for the model.
     * @param is_first If true the first snp is associated, otherwise the second.
     */
    ld_assoc(log_double prior, const arma::vec &alpha, bool is_first);

    /**
     * @see model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

private:
    /**
     * Indicates which snp is associated with the phenotype, if true
     * the first, false the second.
     */
    bool m_is_first;
};

/**
 * This class represents a model where two snps are independently
 * associated to the phenotype.
 *
 * The posterior of this model takes significantly longer than the
 * other models and should only be used on a small subset of
 * snp pairs.
 */
class sindependent
: public model
{
public:
    /**
     * Constructor.
     *
     * @param prior The prior probability for the model.
     * @param num_mc_iterations The number of monte carlo iterations.
     */
    sindependent(log_double prior, const arma::vec &alpha, int num_mc_iterations);

    /**
     * @see model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

private:
    dir_generator m_rdir;

    /**
     * Number of monte carlo iterations.
     */
    int m_num_mc_iterations;
};

#endif /* End of __BAYESIC_MODELS_H__ */
