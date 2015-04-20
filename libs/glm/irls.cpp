#include <iostream>
#include <cfloat>

#include <armadillo>

#include <glm/models/glm_model.hpp>
#include <glm/models/links/glm_link.hpp>
#include <glm/irls.hpp>
#include <dcdflib/libdcdf.hpp>

using namespace arma;

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

    mat Ainv;
    if( pinv( Ainv, A ) )
    {
        return Ainv * ty;
    }
    else
    {
        return vec( );
    }
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
irls(const mat &X, const vec &y, const uvec &missing, const glm_model &model, glm_info &output)
{
    const glm_link &link = model.get_link( );
    vec b = link.init_beta( X, y );
    vec w( X.n_rows );
    vec z( X.n_rows );
    vec eta = X * b;
    vec mu = link.mu( eta );
    vec mu_eta = link.mu_eta( mu );

    int num_iter = 0;
    double old_logl = -DBL_MAX;
    double logl = model.likelihood( mu, y, missing );
    bool invalid_mu = false;
    bool inverse_fail = false;
    while( num_iter < IRLS_MAX_ITERS && ! ( fabs( logl - old_logl ) / ( 0.1 + fabs( logl ) ) < IRLS_TOLERANCE ) )
    {
        w = compute_w( model.var( mu ), mu_eta );
        z = compute_z( eta, mu, mu_eta, y );
        set_missing_to_zero( missing, w );
        b = weighted_least_squares( X, z, w );
        if( b.n_elem <= 0 )
        {
            inverse_fail = true;
            break;
        }
        
        eta = X * b;
        mu = link.mu( eta );
        mu_eta = link.mu_eta( mu );

        if( !model.valid_mu( mu ) )
        {
            invalid_mu = true;
            break;
        }

        old_logl = logl;
        logl = model.likelihood( mu, y, missing );

        num_iter++;
    }

    if( num_iter < IRLS_MAX_ITERS && !invalid_mu && !inverse_fail )
    {
        mat C( b.n_elem, b.n_elem );
        bool inverted = inv( C, X.t( ) *  diagmat( w ) * X );
        if( inverted )
        {
            float dispersion = model.dispersion( mu, y, missing, b.n_elem );
            output.se_beta = sqrt( model.dispersion( mu, y, missing, dispersion ) * diagvec( C ) );
            output.num_iters = num_iter;
            output.converged = true;
            output.success = true;
            output.mu = mu;
            output.logl = model.likelihood( mu, y, missing, dispersion );
            
            vec wald_z = b / output.se_beta;
            vec chi2_value = wald_z % wald_z;
            output.p_value = -1.0 * ones<vec>( chi2_value.n_elem );
            for(int i = 0; i < chi2_value.n_elem; i++)
            {
                try
                {
                    output.p_value[ i ] = 1.0 - chi_square_cdf( chi2_value[ i ], 1 );
                }
                catch(bad_domain_value &e)
                {
                    continue;
                }
            }
        }
        else
        {
            output.success = false;
        }
    }
    else
    {   
        output.num_iters = num_iter;
        output.converged = false;
        output.success = false;
    }

    return b;
}

vec
irls(const mat &X, const vec &y, const glm_model &model, glm_info &output)
{
    uvec missing = zeros<uvec>( y.n_elem );
    return irls( X, y, missing, model, output );
}
