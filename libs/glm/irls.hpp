#ifndef __IRLS_H__
#define __IRLS_H__

#include <armadillo>

#include <models/glm_model.hpp>

/**
 * Maximum number of iterations in the IRLS algorithm.
 */
static const int IRLS_MAX_ITERS = 25;

/**
 * Smallest change in likelihood before terminating the IRLS algorithm.
 */
static double IRLS_TOLERANCE = 10e-8;

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
arma::vec weighted_least_squares(const arma::mat &X, const arma::vec &y, const arma::vec &w);

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
arma::vec irls(const arma::mat &X, const arma::vec &y, const glm_model &model);

#endif /* End of __IRLS_H__ */
