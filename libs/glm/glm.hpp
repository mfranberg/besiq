#ifndef __GLM_H__
#define __GLM_H__

#include <glm/glm_info.hpp>
#include <glm/models/glm_model.hpp>

/**
 * This function fits a generalized linear model. The specific algorithm
 * depends on the model, but in general either matrix inversion for linear
 * regression, and iteratively reweighted least squares for other models.
 *
 * @param X The design matrix (caller is responsible for
 *          adding an intercept).
 * @param y The observations.
 * @param missing Identifies missing sampels by 1 and non-missing by 0.
 * @param model The GLM model to estimate.
 * @param output Output statistics of the estimated betas.
 * @param fast_inversion Use faster but less robust matrix inversion.
 *
 * @return Estimated beta coefficients.
 */
arma::vec glm_fit(const arma::mat &X, const arma::vec &y, const arma::uvec &missing, const glm_model &model, glm_info &output, bool fast_inversion = false);

#endif /* End of __GLM_H__ */
