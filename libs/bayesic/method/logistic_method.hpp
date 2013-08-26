#ifndef __LOGISTIC_METHOD_H__
#define __LOGISTIC_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <models/binomial.hpp>
#include <irls.hpp>
#include <method/method.hpp>
#include <stats/log_scale.hpp>

class logistic_method
: public method_type
{

public:
    logistic_method(method_data_ptr data);
    virtual void run(const snp_row &row1, const snp_row &row2, const std::string &name1, const std::string &name2);

private:
    arma::mat m_design_matrix;
    binomial m_model;
};

#endif /* End of __LOGISTIC_METHOD_H__ */
