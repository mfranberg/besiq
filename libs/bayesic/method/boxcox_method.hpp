#ifndef __BOXCOX_METHOD_H__
#define __BOXCOX_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <glm/glm.hpp>
#include <bayesic/method/method.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/model_matrix.hpp>

/**
 * This class for running the scale invariant method.
 */
class boxcox_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * Lambda goes lambda_start, lambda_start + lambda_step, ..., lambda_end.
     *
     * @param data Additional data required by all methods.
     * @param model_matrix The model matrix.
     * @param is_lm Is this a linear model.
     * @param lambda_start Start value of lambda.
     * @param lambda_end End value of lambda.
     * @param lambda_start Step value of lambda.
     */
    boxcox_method(method_data_ptr data, model_matrix &model_matrix, bool is_lm, float lambda_start, float lambda_end, float lambda_step);
    
    /**
     * Destructor.
     */
    ~boxcox_method();
    
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
     * The included models.
     */
    std::vector<glm_model *> m_model;

    /**
     * The model matrix.
     */
    model_matrix &m_model_matrix;

    /**
     * List of lambda for each link.
     */
    std::vector<float> m_lambda;

    /**
     * Lambda2 to make quantative observations > 0.
     */
    float m_lambda2;
};

#endif /* End of __BOXCOX_METHOD_H__ */
