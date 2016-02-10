#include <glm/models/links/power.hpp>

#include <armadillo>

using namespace arma;

power_link::power_link(float lambda)
    : glm_link::glm_link( "power" )
{
    m_lambda = lambda;
    if( lambda < 1e-6 )
    {
        m_lambda = 0.0;
    }
    else if( lambda > 2.0 )
    {
        m_lambda = 2.0;
    }
}

vec
power_link::init_beta(const mat &X, const vec &y) const
{
    return 0;
}

vec
power_link::mu(const arma::vec &eta) const
{
    if( m_lambda == 0.0 )
    {
        return log( eta );
    }
    else if( m_lambda == 2.0 )
    {
        return exp( eta );
    }
    else if( m_lambda < 1.0 )
    {
        return ( pow( eta, m_lambda ) - 1 ) / m_lambda;
    }
    else
    {
        return pow( 1 + eta*(2-m_lambda), 1/(2-m_lambda) );
    }
}

vec
power_link::eta(const arma::vec &mu) const
{
    if( m_lambda == 0.0 )
    {
        return exp( mu );
    }
    else if( m_lambda == 2.0 )
    {
        return log( mu );
    }
    else if( m_lambda < 1.0 )
    {
        return pow( 1 + mu*m_lambda, 1/m_lambda);
    }
    else
    {
        return (pow( mu, 2 - m_lambda ) - 1) / (2 - m_lambda);
    }
}

vec
power_link::mu_eta(const arma::vec &mu) const
{
    if( m_lambda == 0.0 )
    {
        return exp( mu );
    }
    else if( m_lambda == 2.0 )
    {
        return 1 / mu;
    }
    else if( m_lambda < 1.0 )
    {
        return pow( 1 + mu * m_lambda, 1/m_lambda - 1 );
    }
    else
    {
        return pow( mu, 1 - m_lambda );
    }
}
