#include <models/binomial.hpp>

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
binomial::compute_w(const mat &X, const vec &y, const vec &b) const
{
    vec eta = X * b;
    vec mu = 1.0 / ( 1.0 + exp( -eta ) );
    vec denom = mu % ( 1 - mu );

    vec w = denom % denom % denom;

    return w;
}

vec
binomial::compute_z(const mat &X, const vec &y, const vec &b) const
{
    vec eta = X * b;
    vec mu = 1.0 / ( 1.0 + exp( -eta ) );
    vec z = eta + ( 1.0 / ( mu % ( 1 - mu ) ) ) % ( y - mu );

    return z;
}

double
binomial::likelihood(const mat &X, const vec &y, const vec &b) const
{
    vec eta = X * b;
    vec mu = 1.0 / ( 1.0 + exp( -eta ) );

    return sum( (y % log( mu ) ) + ( (1 - y) % log( 1 - mu ) ) );
}
