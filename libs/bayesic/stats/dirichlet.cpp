#include <assert.h>

#include <stats/dirichlet.hpp>
#include <libdcdf.hpp>

double
dirmult(const arma::vec &x, const arma::vec &alpha)
{
    return exp( ldirmult( x, alpha ) );
}

double
ldirmult(const arma::vec &x, const arma::vec &alpha)
{
    assert( x.n_elem == alpha.n_elem );
    
    double alpha_sum = arma::sum( alpha );
    double x_sum = arma::sum( x );

    double K = lgamma( alpha_sum ) - lgamma( x_sum + alpha_sum );
    double x_alpha = 0.0;

    for(int i = 0; i < x.n_elem; i++)
    {
        x_alpha += lgamma( x[ i ] + alpha[ i ] ) - lgamma( alpha[ i ] );
    }

    return K + x_alpha;
}

double
lbinomial(double n, double k)
{
    assert( k >= 0.0 );
    assert( n >= k );

    if( k == 0.0 )
    {
        return 0.0;
    }

    return lgamma( n + 1.0 ) - lgamma( n - k + 1.0 ) - lgamma( k + 1.0 );
}

dir_generator::dir_generator(unsigned long seed)
    :   m_generator( seed )
{

}

arma::vec
dir_generator::sample(const arma::vec &alpha)
{
    arma::vec x = arma::zeros<arma::vec>( alpha.n_elem );
    for(int i = 0; i < alpha.n_elem; i++)
    {
        double runif =((double) m_generator( )) / ( m_generator.max( ) - m_generator.min( ) );
        x[ i ] = gamma_cdf_inv( runif, alpha[ i ], 1.0 );
    }

    return x / sum( x );
}
