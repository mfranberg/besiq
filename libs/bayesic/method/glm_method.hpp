#ifndef __GLM_METHOD_H__
#define __GLM_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <glm/glm.hpp>
#include <bayesic/method/method.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/model_matrix.hpp>

/**
 * This class is responsible for initializing and repeatedly
 * executing the glm regression on pairs of snps.
 */
class glm_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     */
    glm_method(method_data_ptr data, const glm_model &model, model_matrix &model_matrix);
    
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
     * The glm model used, in this case a binomial model with logit link.
     */
    const glm_model &m_model;

    /**
     * The model matrix that is used.
     */
    model_matrix &m_model_matrix;
};

#endif /* End of __GLM_METHOD_H__ */
