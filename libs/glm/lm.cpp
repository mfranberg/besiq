#include <armadillo>

#include <glm/lm.hpp>
#include <glm/irls.hpp>

using namespace arma;

double loglikelihood(const vec &residuals, double sigma_square, double n)
{
    return -n/2*log(2*datum::pi) - n/2*log( sigma_square ) - 1/(2*sigma_square) * accu( residuals % residuals );
}

vec
lm(const mat &X, const vec &y, const uvec &missing, lm_info &output)
{
    vec w = ones<vec>( y.n_elem );
    set_missing_to_zero( missing, w );
    vec beta = weighted_least_squares( X, y, w );

    double n = accu( w );
    double k = X.n_cols;

    vec mu = X * beta;
    vec residuals = y - mu;
    double sigma_square = as_scalar( trans( residuals ) * ( w % residuals ) / ( n - k ) );
    vec sd = arma::sqrt( sigma_square * diagvec( inv( trans( X ) * ( diagmat( w ) * X ) ) ) );

    output.se_beta = sd;
    output.p_value = chi_square_cdf( beta % beta / ( sd % sd ), 1 );
    output.mu = mu;
    output.logl = loglikelihood( residuals % w, sigma_square, n );

    return beta;
}
