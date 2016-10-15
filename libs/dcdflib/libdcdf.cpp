#include <dcdflib/libdcdf.hpp>
#include <cmath>

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
f_cdf(double x, double d1, double d2)
{
    int which = 1;
    double p;
    double q;
    double x_f = x;
    double d1_f = d1;
    double d2_f = d2;

    int status;
    double bound;

    cdff( &which, &p, &q, &x_f, &d1_f, &d2_f, &status, &bound );

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

double
exp_cdf(double x, double lambda)
{
    if( x < 0 )
    {
        throw bad_domain_value( x );
    }
    if( lambda <= 0 )
    {
        throw bad_domain_value( lambda );
    }

    return 1 - std::exp( -lambda * x );
}
