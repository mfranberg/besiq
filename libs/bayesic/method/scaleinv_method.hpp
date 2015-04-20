#ifndef __SCALEINV_METHOD_H__
#define __SCALEINV_METHOD_H__

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
class scaleinv_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     * @param model_matrix The model matrix.
     */
    scaleinv_method(method_data_ptr data, model_matrix &model_matrix, bool is_lm);
    
    /**
     * Destructor.
     */
    ~scaleinv_method();
    
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
};

#endif /* End of __SCALEINV_METHOD_H__ */
