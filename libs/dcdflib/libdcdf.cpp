#include <stdexcept>

#include <libdcdf.hpp>

extern "C"
{
    #include <cdflib.h>
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
        throw std::runtime_error( "chi_square_cdf: Could not compute chi^2 p-value." );
    }
}
