#include <iostream>
#include <cfloat>

#include <armadillo>

#include <models/glm_model.hpp>
#include <irls.hpp>

using namespace arma;

/**
 * Solves the weighted least square problem:
 *   
 *   W*X*W*b = X*W*y
 *
 * The algorithm uses the singular value decomposition, to
 * compute the solution b:
 *   
 *   b = V * S^-1 * U^t sqrt( w ) * y
 *
 * where wX = U S V^t
 *
 * @param X The design matrix.
 * @param y The right hand side.
 * @param w The weight for each observation.
 *
 * @return The vector b that minimizes the weighted least squares problem.
 */
vec
weighted_least_squares(const mat &X, const vec &y, const vec &w)
{
    /* A = sqrt( w) * X */
    mat A = X;
    A.each_col( ) %= sqrt( w );

    /* ty = sqrt( w ) * y */
    vec ty = y % sqrt( w );

    return pinv( A ) * ty;
}

/**
 * This function performs the iteratively reweighted
 * least squares algorithm to estimate beta coefficients
 * of a genearlized linear model.
 *
 * @param X The design matrix (caller is responsible for
 *          adding an intercept).
 * @param y The observations.
 * @param model The GLM model to estimate.
 *
 * @return Estimated beta coefficients.
 */
vec
irls(const mat &X, const vec &y, const glm_model &model)
{
    vec b = model.init_beta( X, y );
    vec w( X.n_rows );
    vec z( X.n_rows );

    int num_iter = 0;
    double old_logl = -DBL_MAX;
    double logl = model.likelihood( X, y, b );
    while( num_iter < IRLS_MAX_ITERS && ! ( fabs( logl - old_logl ) < IRLS_TOLERANCE ) )
    {
        w = model.compute_w( X, y, b );
        z = model.compute_z( X, y, b );
        b = weighted_least_squares( X, z, w );

        old_logl = logl;
        logl = model.likelihood( X, y, b );
        num_iter++;
    }

    if( num_iter == IRLS_MAX_ITERS + 1 )
    {
        std::cout << "bayesic: warning: IRLS did not converge." << std::endl;
    }

    std::cout << "Converged in " << num_iter << " iterations." << std::endl;

    return b;
}

