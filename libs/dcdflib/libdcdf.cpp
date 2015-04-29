#include <dcdflib/libdcdf.hpp>

extern "C"
{
    #include <dcdflib/cdflib.h>
}

double
chi_square_cdf(double x, unsigned int df)
{
    int which = 1;
    double p;
    double q;
    double x_chi = x;
    double df_chi = df;
    int status;
    double bound;

    cdfchi( &which, &p, &q, &x_chi, &df_chi, &status, &bound );

    if( status == 0 )
    {
        return p;
    }
    else
    {
        throw bad_domain_value( x );
    }
}

double
norm_cdf(double x, double mu, double sd)
{
    int which = 1;
    double p;
    double q;
    double x_norm = x;
    double mu_norm = mu;
    double sd_norm = sd;

    int status;
    double bound;

    cdfnor( &which, &p, &q, &x_norm, &mu_norm, &sd_norm, &status, &bound );

    if( status == 0 )
    {
        return p;
    }
    else
    {
        throw bad_domain_value( x );
    }
}

double
gamma_cdf_inv(double p, double a, double b)
{
    int which = 2;
    double p_gam = p;
    double q = 1.0 - p;
    double x;
    double shape = a;
    double scale = b;
    double bound;
    int status;

    cdfgam( &which, &p_gam, &q, &x, &shape, &scale, &status, &bound );

    if( status == 0 )
    {
        return x;
    }
    else
    {
        throw bad_domain_value( x );
    }
}
