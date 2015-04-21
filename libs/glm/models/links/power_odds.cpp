#include <glm/models/links/power_odds.hpp>

#include <armadillo>

using namespace arma;

power_odds_link::power_odds_link(float lambda)
    : glm_link::glm_link( "power_odds" ),
      m_lambda( lambda )
{
}

vec
power_odds_link::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y+0.5)/2.0;
    vec eta = log( mu / ( 1 - mu ) );
    if( m_lambda != 0.0 )
    {
        eta = ( pow( mu / ( 1 - mu ), m_lambda ) - 1 ) / m_lambda;
    }
    
    return pinv( X ) * eta;
}

vec
power_odds_link::mu(const arma::vec &eta) const
{
    if( m_lambda == 0.0 )
    {
        return exp( eta ) / ( 1 + exp( eta ) );
    }
    else
    {
        return (pow(1 + m_lambda * eta, 1/m_lambda)) / ( 1 + pow(1 + m_lambda * eta, 1/m_lambda));
    }
}

vec
power_odds_link::mu_eta(const arma::vec &mu) const
{
    if( m_lambda == 0.0 )
    {
        return 1.0 / ( mu * ( 1 - mu ) );
    }
    else
    {
        return pow( mu, m_lambda - 1 ) / ( pow( 1-mu, m_lambda + 1 ) );
    }
}
