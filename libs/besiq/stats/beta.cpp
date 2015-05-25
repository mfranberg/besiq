#include <assert.h>

#include <besiq/stats/dirichlet.hpp>

arma::vec
mom_beta(const arma::vec &samples)
{
    double mean = sum( samples ) / samples.n_elem;
    arma::vec dev = samples - mean;

    double var = ( 1.0 / ( samples.n_elem - 1.0 ) ) * sum( dev % dev );

    arma::vec alpha = arma::ones<arma::vec>( 2 );
    if( var < mean * ( 1 - mean ) )
    {
        alpha[ 0 ] = mean * ( ( mean * ( 1 - mean ) / var ) - 1 );
        alpha[ 1 ] = ( 1 - mean ) * ( ( mean * ( 1 - mean ) / var ) - 1 );
    }

    return alpha;
}
