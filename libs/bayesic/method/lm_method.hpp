#ifndef __LINEAR_METHOD_H__
#define __LINEAR_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <glm/lm.hpp>
#include <bayesic/method/method.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/model_matrix.hpp>

/**
 * This class is responsible for initializing and repeatedly
 * executing the linear regression on pairs of snps, this
 * versions treats the SNPs as factor variables instead of
 * ordinals.
 */
class lm_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     * @param model_matrix The model matrix.
     */
    lm_method(method_data_ptr data, model_matrix &model_matrix);
    
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
     * The model matrix.
     */
    model_matrix &m_model_matrix;
};

#endif /* End of __LINEAR_METHOD_H__ */
