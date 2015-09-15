#include <glm/models/links/power.hpp>

#include <armadillo>

using namespace arma;

power_link::power_link(float lambda)
    : glm_link::glm_link( "power" )
{
    m_lambda = lambda;
    if( std::abs( lambda ) < 1e-5 )
    {
        m_lambda = 0.0;
    }
}

vec
power_link::init_beta(const mat &X, const vec &y) const
{
    vec eta;
    if( m_lambda != 0.0 )
    {
        eta = sign( y ) % pow( abs( y ), m_lambda );
    }
    else
    {
        eta = log( y );
    }
    
    return pinv( X ) * eta;
}

vec
power_link::mu(const arma::vec &eta) const
{
    if( m_lambda == 0.0 )
    {
        return exp( eta );
    }
    else
    {
        return sign( eta ) % pow( abs( eta ), ( 1 / m_lambda ) );
    }
}

vec
power_link::eta(const arma::vec &mu) const
{
    if( m_lambda == 0.0 )
    {
        return log( mu );
    }
    else
    {
        return sign( mu ) % pow( abs( mu ), m_lambda );
    }
}

vec
power_link::mu_eta(const arma::vec &mu) const
{
    if( m_lambda == 0.0 )
    {
        return 1.0 / mu;
    }
    else
    {
        return m_lambda * pow( abs( mu ), m_lambda - 1 );
    }
}
