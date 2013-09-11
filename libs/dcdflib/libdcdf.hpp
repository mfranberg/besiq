#ifndef __LIBDCDF_H__
#define __LIBDCDF_H__

/**
 * Computes the probability Pr[ X < x ] where X has a chi^2
 * distribution.
 *
 * @param x Observed chi square value.
 * @param df The degrees of freedom.
 *
 * @return The probability Pr[ X < x ].
 */
double chi_square_cdf(double x, unsigned int df);

/**
 * Computes f^-1(p) where f^-1 is the inverse of the cumulative
 * distribution function for a gamma distributed variable.
 *
 * @param p Is the probability Pr[X <= x] for some x.
 * @param a Shape parameter.
 * @param b Scale parameter.
 *
 * @return The x such that Pr[X <= x] = p.
 */
double gamma_cdf_inv(double p, double a, double b);

#endif /* End of __LIBDCDF_H__ */
