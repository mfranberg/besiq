#include <glm/models/links/odds.hpp>

#include <armadillo>

using namespace arma;

vec
odds_link::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = mu / ( 1.0 - mu );
    
    /* Sometimes the initialization behaves badly for the odds-additive
     * scale by taking a too large step in the first iteration. Always
     * shrink the initial betas some to address this.
    */
    return pinv( X ) * eta / 4.0;
}

vec
odds_link::mu(const arma::vec &eta) const
{
    return eta / ( 1 + eta );
}

vec
odds_link::mu_eta(const arma::vec &mu) const
{
    return 1 / ( (mu - 1) % (mu - 1) );
}
