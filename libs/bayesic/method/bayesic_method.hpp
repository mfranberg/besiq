#ifndef __BAYESIC_METHOD_H__
#define __BAYESIC_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <method/method.hpp>
#include <stats/bayesic_models.hpp>
#include <stats/log_scale.hpp>

class bayesic_method
: public method_type
{

public:
    bayesic_method(method_data_ptr data);
    virtual void run(const snp_row &row1, const snp_row &row2, const std::string &name1, const std::string &name2);

private:
    std::vector<model *> m_models;
    arma::vec m_weight;
};

#endif /* End of __BAYESIC_METHOD_H__ */
