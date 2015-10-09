#ifndef __SEPARATE_METHOD_H__
#define __SEPARATE_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <glm/glm.hpp>
#include <besiq/method/method.hpp>
#include <besiq/stats/log_scale.hpp>
#include <besiq/model_matrix.hpp>

/**
 * This class for running each interaction parameter separately.
 */
class separate_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     */
    separate_method(method_data_ptr data, glm_model *model);
    
    /**
     * Destructor.
     */
    ~separate_method();
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();

    /**
     * @see method_type::run.
     */
    virtual double run(const snp_row &row1, const snp_row &row2, float *output);

private:
    /**
     * The included models.
     */
    std::vector<model_matrix *> m_model_matrix;

    /**
     * The glm model.
     */
    glm_model *m_model;

};

#endif /* End of __SEPARATE_METHOD_H__ */
