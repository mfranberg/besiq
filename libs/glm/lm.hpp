#ifndef __LM_H__
#define __LM_H__

#include <armadillo>

#include <glm/glm_info.hpp>
#include <glm/models/glm_model.hpp>

/**
 * Computes the log-likelihood of the linear model.
 *
 * @param residuals The residual vector, zero if missing.
 * @param sigma_square Estimate of the variance.
 * @param n The number of non-missing individuals.
 *
 * @return The log-likelihood of the linear model.
 */
double loglikelihood(const arma::vec &residuals, double sigma_square, double n);

/**
 * This function solves the linear least squares problem.
 *
 * @param X The design matrix (caller is responsible for
 *          adding an intercept).
 * @param y The observations.
 * @param missing Identifies missing sampels by 1 and non-missing by 0.
 * @param output Output statistics of the estimated betas.
 *
 * @return Estimated beta coefficients.
 */
arma::vec lm(const arma::mat &X, const arma::vec &y, const arma::uvec &missing, const glm_model &model, glm_info &output);

#endif /* End of __LM_H__ */
