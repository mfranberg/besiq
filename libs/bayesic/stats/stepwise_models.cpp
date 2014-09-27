#include <bayesic/stats/stepwise_models.hpp>

using namespace arma;

lm_full::lm_full()
: stepwise_model::stepwise_model( 0 )
{

}

log_double
lm_full::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    arma::mat count = joint_count_cont( row1, row2, phenotype, weight );
    arma::vec mu_full = count.col( 0 ) / count.col( 1 );
    arma::vec residual = count.col( 1 ) % mu_full % mu_full - 2 * mu_full % count.col( 0 ) + count.col( 2 );
    double k = 9;
    double n = sum( count.col( 1 ) );
    double sigma_square = sum( residual ) / ( n - k );

    *is_valid = arma::min( count.col( 1 ) ) >= 5;
    return log_double::from_log( -(n/2)*log(2*datum::pi) - (n/2)*log( sigma_square ) - 1/(2*sigma_square) * sum( residual ) );
}

intercept::intercept()
: stepwise_model::stepwise_model( 8 )
{

}

log_double
intercept::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    arma::mat count = joint_count_cont( row1, row2, phenotype, weight );

    double mu = sum( count.col( 0 ) ) / sum( count.col( 1 ) );
    double residual = sum( count.col( 1 ) ) * mu * mu - 2 * mu * sum( count.col( 0 ) ) + sum( count.col( 2 ) );
    double k = 1;
    double n = sum( count.col( 1 ) );
    double sigma_square = residual / ( n - k );

    *is_valid = arma::min( count.col( 1 ) ) >= 5;
    return log_double::from_log( -(n/2)*log(2*datum::pi) - (n/2)*log( sigma_square ) - 1/(2*sigma_square) * residual );
}

single::single(bool is_first)
: stepwise_model::stepwise_model( 6 ),
  m_is_first( is_first )
{

}

log_double
single::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    const snp_row &snp1 = m_is_first ? row1 : row2;
    const snp_row &snp2 = m_is_first ? row2 : row1;
    
    arma::mat count = joint_count_cont( snp1, snp2, phenotype, weight );
    arma::mat snp_pheno = zeros<mat>( 3, 3 );
    for(int i = 0; i < 3; i++)
    {
        snp_pheno( i, 0 ) = count( 3*i, 0 ) + count( 3*i + 1, 0 ) + count( 3*i + 2, 0 );
        snp_pheno( i, 1 ) = count( 3*i, 1 ) + count( 3*i + 1, 1 ) + count( 3*i + 2, 1 );
        snp_pheno( i, 2 ) = count( 3*i, 2 ) + count( 3*i + 1, 2 ) + count( 3*i + 2, 2 );
    }
    
    arma::vec mu_single = snp_pheno.col( 0 ) / snp_pheno.col( 1 );
    arma::vec residual = snp_pheno.col( 1 ) % mu_single % mu_single - 2 * mu_single % snp_pheno.col( 0 ) + snp_pheno.col( 2 );
    double k = 3;
    double n = sum( snp_pheno.col( 1 ) );
    double sigma_square = sum( residual ) / ( n - k );

    *is_valid = arma::min( count.col( 1 ) ) >= 5;
    return log_double::from_log( -(n/2)*log(2*datum::pi) - (n/2)*log( sigma_square ) - 1/(2*sigma_square) * sum( residual ) );
}
