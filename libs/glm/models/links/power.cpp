#include <glm/models/links/power.hpp>

#include <armadillo>

using namespace arma;

power_link::power_link(float lambda, float lambda2)
    : glm_link::glm_link( "power" ),
      m_lambda( lambda ),
      m_lambda2( lambda2 )
{
}

vec
power_link::init_beta(const mat &X, const vec &y) const
{
    vec mu = y;
    vec eta = log( mu + m_lambda2 );

    if( m_lambda != 0.0 )
    {
        eta = pow( mu + m_lambda2, m_lambda );
    }
    
    return pinv( X ) * eta;
}

vec
power_link::mu(const arma::vec &eta) const
{
    if( m_lambda == 0.0 )
    {
        return exp( eta ) - m_lambda2;
    }
    else
    {
        return pow( eta, ( 1 / m_lambda ) ) - m_lambda2;
    }
}

vec
power_link::mu_eta(const arma::vec &mu) const
{
    if( m_lambda == 0.0 )
    {
        return 1.0 / ( mu + m_lambda2 );
    }
    else
    {
        return m_lambda * pow( mu + m_lambda2, ( m_lambda - 1 ) );
    }
}
