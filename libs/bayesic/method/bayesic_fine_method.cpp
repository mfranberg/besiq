#include <bayesic/method/bayesic_fine_method.hpp>

bayesic_fine_method::bayesic_fine_method(method_data_ptr data, int num_mc_iterations, arma::vec alpha)
    : bayesic_method::bayesic_method( data )
{
    double int_prior = 1.0 / ( 2.0 * data->num_interactions );
    
    if( data->single_prior > 0.0 )
    {
        int_prior = data->single_prior;
    }

    std::vector<model *> models;
    models.push_back( new saturated( int_prior, alpha ) );
    models.push_back( new sindependent( 1.0 - int_prior, alpha, num_mc_iterations ) );

    bayesic_method::set_models( models );
}
