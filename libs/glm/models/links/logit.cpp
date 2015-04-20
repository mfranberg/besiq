#include <glm/models/links/logit.hpp>

#include <armadillo>

using namespace arma;

logit_link::logit_link()
    : glm_link::glm_link( "logit" )
{
}

vec
logit_link::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = log( mu / ( 1.0 - mu ) );
    
    return pinv( X ) * eta;
}

vec
logit_link::mu(const arma::vec &eta) const
{
    return 1.0 / ( 1.0 + exp( -eta ) );
}

vec
logit_link::mu_eta(const arma::vec &mu) const
{
    return 1.0 / ( mu % ( 1.0 - mu ) );
}
