#include <glm/models/logcomplement.hpp>

#include <armadillo>

using namespace arma;

vec
logcomplement::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = log( 1 - mu );
    
    return pinv( X ) * eta;
}

vec
logcomplement::mu(const arma::vec &eta) const
{
    return 1.0 - exp( eta );
}

vec
logcomplement::mu_eta(const arma::vec &mu) const
{
    return -1 / ( 1.0 - mu );
}

vec
logcomplement::var(const arma::vec &mu) const
{
    return mu % ( 1.0 - mu );
}

double
logcomplement::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const
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
