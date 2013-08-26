#include <assert.h>

#include <stats/dirichlet.hpp>

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
