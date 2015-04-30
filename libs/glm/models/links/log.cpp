#include <glm/models/links/log.hpp>

#include <armadillo>

using namespace arma;

log_link::log_link()
    : glm_link::glm_link( "log" )
{
}

vec
log_link::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = log( mu );
    
    return pinv( X ) * eta;
}

vec
log_link::mu(const arma::vec &eta) const
{
    return exp( eta );
}

vec
log_link::eta(const arma::vec &mu) const
{
    return log( mu );
}

vec
log_link::mu_eta(const arma::vec &mu) const
{
    return 1.0 / mu;
}
