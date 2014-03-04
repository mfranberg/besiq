#include <bayesic/stats/loglinear_models.hpp>

using namespace arma;

full::full()
: loglinear_model::loglinear_model( 16, 0 )
{

}

log_double
full::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    arma::mat count = joint_count( row1, row2, phenotype, weight ) + 1.0;
    arma::mat p_full = count / accu( count );
    
    *is_valid = arma::min( arma::min( count ) ) >= 5;
    return log_double::from_log( accu( count % arma::log( p_full ) ) );
}

block::block()
: loglinear_model::loglinear_model( 8, 8 )
{

}

log_double
block::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    arma::mat counts = joint_count( row1, row2, phenotype, weight ) + 1.0;

    arma::vec snp_snp = sum( counts, 1 );
    arma::rowvec pheno = sum( counts, 0 );

    double n = accu( counts );

    double likelihood = 0.0;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            for(int k = 0; k < 2; k++)
            {
                double n_ijk = counts( 3 * i + j, k );
                double p = ( snp_snp( 3 * i + j ) * pheno( k ) / ( n * n ) );
                likelihood += n_ijk * log( p );
            }
        }
    }

    *is_valid = arma::min( snp_snp ) >= 5 && arma::min( pheno ) >= 5;
    return log_double::from_log( likelihood );
}

partial::partial(bool is_first)
: loglinear_model::loglinear_model( 10, 6 ),
  m_is_first( is_first )
{

}

log_double
partial::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    const snp_row &snp1 = m_is_first ? row1 : row2;
    const snp_row &snp2 = m_is_first ? row2 : row1;
    
    arma::mat counts = joint_count( snp1, snp2, phenotype, weight ) + 1.0;
    double n = accu( counts );
    arma::vec snp_snp = sum( counts, 1 );
    arma::mat snp_pheno = zeros<mat>( 3, 2 );
    for(int i = 0; i < 3; i++)
    {
        snp_pheno( i, 0 ) = counts( 3*i, 0 ) + counts( 3*i + 1, 0 ) + counts( 3*i + 2, 0 );
        snp_pheno( i, 1 ) = counts( 3*i, 1 ) + counts( 3*i + 1, 1 ) + counts( 3*i + 2, 1 );
    }
    
    double likelihood = 0.0;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            for(int k = 0; k < 2; k++)
            {
                double n_ijk = counts( 3 * i + j, k );
                double p = ( snp_snp[ 3 * i + j ] * snp_pheno( i, k ) ) / ( ( snp_pheno( i, 0 ) + snp_pheno( i, 1 ) ) * n );
                likelihood += n_ijk * log( p );
            }
        }
    }

    *is_valid = arma::min( snp_snp ) >= 5 && arma::min( arma::min( snp_pheno ) ) >= 5;
    return log_double::from_log( likelihood );
}
