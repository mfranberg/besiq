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

#endif /* End of __LIBDCDF_H__ */
