#include <glm/models/binomial.hpp>

#include <armadillo>

using namespace arma;

binomial::binomial(const std::string &link_name)
    : glm_model( "binomial", link_name )
{
}

binomial::binomial(glm_link *link)
    : glm_model( "binomial", link )
{
}

binomial::~binomial()
{
}

bool
binomial::valid_mu(const arma::vec &mu) const
{
    for(int i = 0; i < mu.n_elem; i++)
    {
        if( mu[ i ] <= 0.0 || mu[ i ] >= 1.0 )
        {
            return false;
        }
    }
    
    return true;
}

vec
binomial::var(const arma::vec &mu) const
{
    return mu % ( 1.0 - mu );
}

double
binomial::dispersion(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float k) const
{
    return 1.0;
}

double
binomial::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing, float dispersion) const
{
    double loglikelihood = 0.0;
    for(int i = 0; i < y.n_elem; i++)
    {
        if( missing[ i ] == 0 )
        {
            loglikelihood += y[ i ] * log( mu[ i ] ) + ( 1 - y[ i ] ) * log( 1 - mu[ i ] );
        }
    }

    return loglikelihood;
}

bool
binomial::is_binary() const
{
    return true;
}
