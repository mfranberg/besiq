#ifndef __LIBDCDF_H__
#define __LIBDCDF_H__

#include <sstream>
#include <stdexcept>

class bad_domain_value: public std::exception
{
public:
    /**
     * Constructor.
     *
     * @param value The value that caused an error when computing a p-value.
     */
    bad_domain_value(double value)
    {
        std::ostringstream error_message;
        error_message << "Could not compute p-value for: " << value << std::endl;
        m_message = error_message.str( ).c_str( );
    }

    /**
     * Destructor.
     */
    virtual ~bad_domain_value() throw()
    {
    }

    /**
     * Returns the error message.
     *
     * @return the error message.
     */
    virtual const char* what() const throw()
    {
        return m_message.c_str( );
    }

private:
    /**
     * The error message.
     */
    std::string m_message;
};

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
 * Computes the probability Pr[ X < x ] where X has a normal
 * distribution.
 *
 * @param x Observed normal value.
 * @param mu The mean value.
 * @param sd The standard deviation.
 *
 * @return The probability Pr[ X < x ].
 */
double norm_cdf(double x, double mu, double sd);

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
