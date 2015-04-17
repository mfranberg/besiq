#include <bayesic/stats/normal_models.hpp>

using namespace arma;

normal_full::normal_full()
: closed_form_model::closed_form_model( 9 )
{

}

log_double
normal_full::prob(const arma::mat &count)
{
    arma::vec mu_full = count.col( 0 ) / count.col( 1 );
    arma::vec residual = count.col( 1 ) % mu_full % mu_full - 2 * mu_full % count.col( 0 ) + count.col( 2 );
    double k = 9;
    double n = sum( count.col( 1 ) );
    double sigma_square = sum( residual ) / ( n - k );

    return log_double::from_log( -(n/2)*log(2*datum::pi) - (n/2)*log( sigma_square ) - 1/(2*sigma_square) * sum( residual ) );
}

normal_null::normal_null()
: closed_form_model::closed_form_model( 1 )
{

}

log_double
normal_null::prob(const arma::mat &count)
{
    double mu = sum( count.col( 0 ) ) / sum( count.col( 1 ) );
    double residual = sum( count.col( 1 ) ) * mu * mu - 2 * mu * sum( count.col( 0 ) ) + sum( count.col( 2 ) );
    double k = 1;
    double n = sum( count.col( 1 ) );
    double sigma_square = residual / ( n - k );

    return log_double::from_log( -(n/2)*log(2*datum::pi) - (n/2)*log( sigma_square ) - 1/(2*sigma_square) * residual );
}

normal_single::normal_single(bool is_first)
: closed_form_model::closed_form_model( 3 ),
  m_is_first( is_first )
{

}

log_double
normal_single::prob(const arma::mat &count)
{
    arma::mat snp_pheno = zeros<mat>( 3, 3 );
    for(int i = 0; i < 3; i++)
    {
        if( m_is_first )
        {
            snp_pheno( i, 0 ) = count( 3*i, 0 ) + count( 3*i + 1, 0 ) + count( 3*i + 2, 0 );
            snp_pheno( i, 1 ) = count( 3*i, 1 ) + count( 3*i + 1, 1 ) + count( 3*i + 2, 1 );
            snp_pheno( i, 2 ) = count( 3*i, 2 ) + count( 3*i + 1, 2 ) + count( 3*i + 2, 2 );
        }
        else
        {
            snp_pheno( i, 0 ) = count( 3*0 + i, 0 ) + count( 3*1 + i, 0 ) + count( 3*2 + i, 0 );
            snp_pheno( i, 1 ) = count( 3*0 + i, 1 ) + count( 3*1 + i, 1 ) + count( 3*2 + i, 1 );
            snp_pheno( i, 2 ) = count( 3*0 + i, 2 ) + count( 3*1 + i, 2 ) + count( 3*2 + i, 2 );
        }
    }
    
    arma::vec mu_single = snp_pheno.col( 0 ) / snp_pheno.col( 1 );
    arma::vec residual = snp_pheno.col( 1 ) % mu_single % mu_single - 2 * mu_single % snp_pheno.col( 0 ) + snp_pheno.col( 2 );
    double k = 3;
    double n = sum( snp_pheno.col( 1 ) );
    double sigma_square = sum( residual ) / ( n - k );

    return log_double::from_log( -(n/2)*log(2*datum::pi) - (n/2)*log( sigma_square ) - 1/(2*sigma_square) * sum( residual ) );
}
