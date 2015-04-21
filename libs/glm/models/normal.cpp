#include <glm/models/normal.hpp>

#include <armadillo>

using namespace arma;

normal::normal(const std::string &link_name)
    : glm_model( "normal", link_name )
{
}

normal::normal(glm_link *link)
    : glm_model( "normal", link )
{
}

normal::~normal()
{
}

bool
normal::valid_mu(const arma::vec &mu) const
{ 
    return mu.is_finite( );
}

vec
normal::var(const arma::vec &mu) const
{
    return arma::ones<arma::vec>( mu.n_elem );
}

double 
normal::dispersion(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float k) const
{
    return arma::as_scalar( ( arma::trans( (y - mu) % ( y - mu ) ) * ( 1 - missing ) ) / ( arma::sum( 1 - missing ) - k ) );
}

double
normal::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float dispersion) const
{
    float sigma2 = dispersion;
    return arma::as_scalar( ( -0.5*log( 2 * datum::pi ) -0.5*log( sigma2 ) - (1.0/(2*sigma2))*arma::trans( ( ( y - mu ) % ( y - mu ) ) ) ) * ( 1 - missing ) );
}
