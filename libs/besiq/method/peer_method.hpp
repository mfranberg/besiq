#ifndef __PEER_METHOD_H__
#define __PEER_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <besiq/method/method.hpp>
#include <besiq/method/wald_method.hpp>
#include <besiq/stats/log_scale.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing the case-only test by Lewinger et al.
 */
class peer_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    peer_method(method_data_ptr data);
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();

    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, float *output);

private: 
    void compute_ld_p(const arma::mat &counts, float *ld_case_z, float *ld_contrast_z);
    arma::mat encode_counts(const arma::mat &counts, size_t snp1_onestart, size_t snp2_onestart);
    /**
     * A weight > 0 associated with each sample, that allows for
     * covariate adjustment.
     */
    arma::vec m_weight;
};

#endif /* End of __PEER_METHOD_H__ */
