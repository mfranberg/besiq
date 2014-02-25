#include <glm/models/binomial.hpp>

#include <armadillo>

using namespace arma;

vec
binomial::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = mu / ( 1.0 - mu );
    
    return pinv( X ) * eta;
}

vec
binomial::mu(const arma::vec &eta) const
{
    return 1.0 / ( 1.0 + exp( -eta ) );
}

bool
binomial::valid_mu(const arma::vec &mu) const
{
    for(int i = 0; i < mu.n_elem; i++)
    {
        if( mu[ i ] < 0.0 || mu[ i ] > 1.0 )
        {
            return false;
        }
    }
    
    return true;
}

vec
binomial::mu_eta(const arma::vec &mu) const
{
    return 1.0 / ( mu % ( 1.0 - mu ) );
}

vec
binomial::var(const arma::vec &mu) const
{
    return mu % ( 1.0 - mu );
}

double
binomial::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const
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
    //return sum( (y % log( mu ) ) + ( (1 - y) % log( 1 - mu ) ) );
}
