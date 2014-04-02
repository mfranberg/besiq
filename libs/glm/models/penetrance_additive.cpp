#include <glm/models/penetrance_additive.hpp>

#include <armadillo>

using namespace arma;

vec
penetrance_additive::init_beta(const mat &X, const vec &y) const
{
    vec mu = (y + 0.5) / 2.0;
    vec eta = mu / ( 1.0 - mu );
    
    return 0.1 * arma::ones<arma::vec>( X.n_cols );
}

vec
penetrance_additive::mu(const arma::vec &eta) const
{
    return eta;
}

bool
penetrance_additive::valid_mu(const arma::vec &mu) const
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
penetrance_additive::mu_eta(const arma::vec &mu) const
{
    return arma::ones<arma::vec>( mu.n_elem );
}

vec
penetrance_additive::var(const arma::vec &mu) const
{
    return mu % ( 1.0 - mu );
}

double
penetrance_additive::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const
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
