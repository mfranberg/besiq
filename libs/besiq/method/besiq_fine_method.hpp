#ifndef __BAYESIC_FINE_METHOD_H__
#define __BAYESIC_FINE_METHOD_H__

#include <string>
#include <vector>

#include <besiq/method/besiq_method.hpp>
#include <besiq/stats/log_scale.hpp>
#include <plink/snp_row.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing the besiq method on pairs of snps.
 */
class besiq_fine_method
: public besiq_method
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     * @param num_mc_iterations The number of monte carlo iterations.
     */
    besiq_fine_method(method_data_ptr data, int num_mc_iterations, arma::vec alpha = arma::ones<arma::vec>( 2 ));
};

#endif /* End of __BAYESIC_FINE_METHOD_H__ */
