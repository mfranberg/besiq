#ifndef __NORMAL_MODELS_H__
#define __NORMAL_MODELS_H__

#include <armadillo>

#include <besiq/stats/log_scale.hpp>
#include <besiq/stats/closed_form_models.hpp>

/**
 * This model represents a full model in the sense
 * that each cell in the table has it is own parameter.
 */
class normal_full
: public closed_form_model
{
public:
    /**
     * Constructor.
     */
    normal_full();

    /**
     * @see closed_form_model::prob.
     */
    virtual log_double prob(const arma::mat &count);
};

/**
 * This class represents a model where no snp is associated with
 * the phenotype. The snps may however be in ld.
 */
class normal_null
: public closed_form_model
{
public:
    /**
     * Constructor.
     */
    normal_null();

    /**
     * @see closed_form_model::prob.
     */
    virtual log_double prob(const arma::mat &count);
};

/**
 * This class represents a model where one snp is associated with
 * the phenotype, and the other possibly in ld with the first.
 */
class normal_single
: public closed_form_model
{
public:
    /**
     * Constructor.
     *
     * @param is_first If true the first snp is associated, otherwise the second.
     */
    normal_single(bool is_first);

    /**
     * @see closed_form_model::prob.
     */
    virtual log_double prob(const arma::mat &count);

private:
    /**
     * Indicates which snp is associated with the phenotype, if true
     * the first, false the second.
     */
    bool m_is_first;
};

#endif /* End of __NORMAL_MODELS_H__ */
