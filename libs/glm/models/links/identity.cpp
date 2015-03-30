#include <glm/models/links/identity.hpp>

#include <armadillo>

using namespace arma;

vec
identity_link::init_beta(const mat &X, const vec &y) const
{ 
    vec mu = (y + 0.5) / 2.0;
    vec eta = mu;
    
    return pinv( X ) * eta;
}

vec
identity_link::mu(const arma::vec &eta) const
{
    return eta;
}

vec
identity_link::mu_eta(const arma::vec &mu) const
{
    return arma::ones<arma::vec>( mu.n_elem );
}
