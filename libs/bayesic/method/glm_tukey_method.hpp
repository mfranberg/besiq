#ifndef __LOGISTIC_TUKEY_METHOD_H__
#define __LOGISTIC_TUKEY_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <glm/models/glm_model.hpp>
#include <glm/irls.hpp>
#include <bayesic/method/method.hpp>
#include <bayesic/stats/log_scale.hpp>

/**
 * This class is responsible for initializing and repeatedly
 * executing the glm regression on pairs of snps, this
 * versions treats the SNPs as factor variables instead of
 * ordinals, but have a single parameter for the interaction.
 */
class glm_tukey_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     */
    glm_tukey_method(method_data_ptr data, const glm_model &model);
    
    /**
     * @see method_type::init.
     */
    virtual void init(std::ostream &output);

    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, std::ostream &output);

private:
    /**
     * The alternative design matrix, interactions effects are present.
     *
     * The first 5 columns are:
     * snp1_1, snp1_2, snp2_1, snp2_2, snp1_1 * snp2_1 or snp1_1 * snp2_2 or snp1_2 * snp2_1 or
     * snp1_2 * snp2_2 and intercept. The rest are covariates.
     */
    arma::mat m_alt_design_matrix;

    /**
     * The null design matrix, only additive effects are present.
     *
     * The first 5 columns are:
     * snp1_1, snp1_2, snp2_1, snp2_2 and intercept. The rest are covariates.
     */
    arma::mat m_null_design_matrix;

    /**
     * The glm model used, in this case a binomial model with logit link.
     */
    const glm_model &m_model;

    /**
     * Updates the design matrix for the alternative (saturated) model.
     *
     * @param row1 First snp.
     * @param row2 Second snp.
     * @param design_matrix The design matrix.
     * @param missing The missing samples.
     */
    void init_alt(const snp_row &row1, const snp_row &row2, arma::mat &design_matrix, arma::uvec &missing);
    
    /**
     * Updates the design matrix for the null model.
     *
     * @param row1 First snp.
     * @param row2 Second snp.
     * @param design_matrix The design matrix.
     * @param missing The missing samples.
     */
    void init_null(const snp_row &row1, const snp_row &row2, arma::mat &design_matrix, arma::uvec &missing);
};

#endif /* End of __LOGISTIC_TUKEY_METHOD_H__ */
