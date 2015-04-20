#ifndef __GLM_INFO_H__
#define __GLM_INFO_H__

#include <armadillo>

/**
 * Information about the estimated parameters.
 */
struct glm_info
{
    /**
    * Standard error of estimated beta.
    */
    arma::vec se_beta;

    /**
    * P-value for each beta, based on Wald test.
    */
    arma::vec p_value;
    
    /**
     * Estimated mean value for each individual.
     */
    arma::vec mu;
    
    /**
     * Log likelihood of the model.
     */
    double logl;
    
    /**
     * Converged and had no inversion problems.
     */
    bool success;

    /**
    * For iterative algorithms, the number of iterations.
    */
    unsigned int num_iters;
    
    /**
    * For iterative algorithms, True if converged, false otherwise.
    */
    bool converged;
};

#endif /* End of __GLM_INFO_H__ */
