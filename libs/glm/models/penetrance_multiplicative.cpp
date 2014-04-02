#include <glm/models/penetrance_multiplicative.hpp>

#include <armadillo>

using namespace arma;

vec
penetrance_multiplicative::init_beta(const mat &X, const vec &y) const
{
    
    return -0.4 * arma::ones<arma::vec>( X.n_cols );
}

vec
penetrance_multiplicative::mu(const arma::vec &eta) const
{
    return exp( eta );
}

bool
penetrance_multiplicative::valid_mu(const arma::vec &mu) const
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
penetrance_multiplicative::mu_eta(const arma::vec &mu) const
{
    return 1.0 / mu;
}

vec
penetrance_multiplicative::var(const arma::vec &mu) const
{
    return mu % ( 1.0 - mu );
}

double
penetrance_multiplicative::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const
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
