#include <glm/models/normal.hpp>

#include <armadillo>

using namespace arma;

normal::normal(const std::string &link_name)
    : glm_model( "normal", link_name )
{
}

normal::~normal()
{
}

bool
normal::valid_mu(const arma::vec &mu) const
{ 
    return true;
}

vec
normal::var(const arma::vec &mu) const
{
    return arma::ones<arma::vec>( mu.n_elem );
}

double 
normal::dispersion(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float k) const
{
    return as_scalar( ( ( (y - mu) % ( y - mu ) ) % ( 1 - missing ) ) / ( sum( missing ) - k ) );
}

double
normal::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const
{
    return as_scalar( (-1.0/2*log(2*datum::pi) - 1.0/2*log( 1.0 ) - 1/(2*1.0) * ((y - mu) % (y - mu ))) * (1 - missing) );
}
