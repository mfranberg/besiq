#ifndef __BAYES_FAST_METHOD_H__
#define __BAYES_FAST_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <besiq/method/method.hpp>
#include <besiq/stats/besiq_models.hpp>
#include <besiq/stats/log_scale.hpp>
#include <plink/snp_row.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing the besiq method on pairs of snps.
 */
class bayes_fast_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    bayes_fast_method(method_data_ptr data, arma::vec alpha = arma::ones<arma::vec>( 2 ));
    virtual ~bayes_fast_method();

    /**
     * Sets the models that should be compared. Memory will be
     * managed by this class.
     *
     * @param models A list of models.
     */
    void set_models(const std::vector<model *> &models);

    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();
    
    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, float *output);

private:
    /**
     * The different models that is part of the besiq method,
     * saturated, ld and null.
     */
    std::vector<model *> m_models;
    
    /**
     * A weight > 0 associated with each sample, that allows for
     * covariate adjustment.
     */
    arma::vec m_weight;
};

#endif /* End of __BAYES_FAST_METHOD_H__ */
