#include <bayesic/method/bayesic_fine_method.hpp>

bayesic_fine_method::bayesic_fine_method(method_data_ptr data, int num_mc_iterations, arma::vec alpha)
    : bayesic_method::bayesic_method( data )
{
    std::vector<model *> models;
    models.push_back( new saturated( 1.0 / 2.0, alpha ) );
    models.push_back( new sindependent( 1.0 / 2.0, alpha, num_mc_iterations ) );


    bayesic_method::set_models( models );
}
