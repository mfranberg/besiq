#include <iostream>
#include <cfloat>

#include <armadillo>

#include <models/glm_model.hpp>
#include <irls.hpp>
#include <libdcdf.hpp>

using namespace arma;

/**
* Sets the weights of missing observations to zero, so that
* the will not influence the regression.
*/
void
set_missing_to_zero(const uvec &missing, vec &w)
{   
    for(int i = 0; i < w.size( ); i++)
    {
        if( missing[ i ] == 1 )
        {
            w[ i ] = 0.0;
        }
    }
}

vec
chi_square_cdf(const vec &x, unsigned int df)
{
    vec p = ones<vec>( x.n_elem );
    for(int i = 0; i < x.n_elem; i++)
    {
        p[ i ] = chi_square_cdf( x[ i ], df );
    }

    return p;
}

vec
weighted_least_squares(const mat &X, const vec &y, const vec &w)
{
    /* A = sqrt( w) * X */
    mat A = diagmat( sqrt( w ) ) * X;

    /* ty = sqrt( w ) * y */
    vec ty = y % sqrt( w );

    return pinv( A ) * ty;
}

vec
compute_z(const vec &eta, const vec &mu, const vec &mu_eta, const vec &y)
{
    return eta + mu_eta % ( y - mu );
}

vec
compute_w(const vec &var, const vec& mu_eta)
{
    return 1.0 / ( var % ( mu_eta % mu_eta ) );
}

vec
irls(const mat &X, const vec &y, const uvec &missing, const glm_model &model, irls_info &output)
{
    vec b = model.init_beta( X, y );
    vec w( X.n_rows );
    vec z( X.n_rows );
    vec eta = X * b;
    vec mu = model.mu( eta );
    vec mu_eta = model.mu_eta( mu );

    int num_iter = 0;
    double old_logl = -DBL_MAX;
    double logl = model.likelihood( mu, y, missing );
    while( num_iter < IRLS_MAX_ITERS && ! ( fabs( logl - old_logl ) < IRLS_TOLERANCE ) )
    {
        w = compute_w( model.var( mu ), mu_eta );
        z = compute_z( eta, mu, mu_eta, y );
        set_missing_to_zero( missing, w );
        b = weighted_least_squares( X, z, w );
        
        eta = X * b;
        mu = model.mu( eta );
        mu_eta = model.mu_eta( mu );

        old_logl = logl;
        logl = model.likelihood( mu, y, missing );

        num_iter++;
    }

    if( num_iter <= IRLS_MAX_ITERS )
    {
        mat C = inv( X.t( ) *  diagmat( w ) * X );
        output.se_beta = sqrt( diagvec( C ) );
        vec wald_z = b / output.se_beta;
        output.p_value = 1.0 - chi_square_cdf( wald_z % wald_z, 1 );
        output.num_iters = num_iter;
        output.converged = true;
        output.mu = mu;
    }
    else
    {   
        output.num_iters = num_iter;
        output.converged = false;
    }

    return b;
}

vec
irls(const mat &X, const vec &y, const glm_model &model, irls_info &output)
{
    uvec missing = zeros<uvec>( y.n_elem );
    return irls( X, y, missing, model, output );
}
