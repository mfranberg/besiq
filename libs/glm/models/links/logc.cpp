#include <glm/models/links/logc.hpp>

#include <armadillo>

using namespace arma;

logc_link::logc_link()
    : glm_link::glm_link( "logc" )
{
}

vec
logc_link::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = log( 1 - mu );
    
    return pinv( X ) * eta;
}

vec
logc_link::mu(const arma::vec &eta) const
{
    return 1.0 - exp( eta );
}

vec
logc_link::mu_eta(const arma::vec &mu) const
{
    return -1 / ( 1.0 - mu );
}
