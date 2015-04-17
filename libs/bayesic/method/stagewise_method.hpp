#ifndef __STAGEWISE_METHOD_H__
#define __STAGEWISE_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <bayesic/method/method.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/stats/closed_form_models.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing a stepwise method where p-values for a series
 * of more complex models are computed.
 */
class stagewise_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    stagewise_method(method_data_ptr data, const std::string &model);
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();

    /**
     * Returns the number of usable samples.
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     * @param phenotype The phenotype.
     * 
     * @return The number of usable samples.
     */    
    unsigned int num_ok_samples(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype);

    /**
     * Computes a likelihood ratio statistics that measures
     * whether two SNPs are in linkage disequilibrium or not. It
     * should follow a chi^2(1) distribution if there is no LD.
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     *
     * @return The estimated likelihood ratio statistic.
     */
    double compute_ld_lr(const snp_row &row1, const snp_row &row2);

    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, float *output);

private:
    /**
     * Type of model.
     */
    std::string m_model;
    /**
     * A weight > 0 associated with each sample, that allows for
     * covariate adjustment.
     */
    arma::vec m_weight;

    /**
     * The models used.
     */
    std::vector<closed_form_model *> m_models;

    /**
     * Number of usable samples in the last call to run.
     */
    unsigned int m_num_usable_samples;
};

#endif /* End of __STAGEWISE_METHOD_H__ */
