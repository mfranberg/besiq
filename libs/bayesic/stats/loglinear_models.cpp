#include <bayesic/stats/loglinear_models.hpp>

using namespace arma;

full::full()
: loglinear_model::loglinear_model( 16, 0 )
{

}

log_double
full::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    arma::mat count = joint_count( row1, row2, phenotype, weight ) + 0.5;
    arma::vec p_full = count.col( 1 ) / ( count.col( 0 ) + count.col( 1 ) );

    *is_valid = arma::min( arma::min( count ) ) >= 5;
    return log_double::from_log( accu( count.col( 1 ) % arma::log( p_full ) + count.col( 0 ) % arma::log( 1 - p_full ) ) );
}

block::block()
: loglinear_model::loglinear_model( 8, 8 )
{

}

log_double
block::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid)
{
    arma::mat counts = joint_count( row1, row2, phenotype, weight ) + 0.5;

    arma::rowvec pheno = sum( counts, 0 );
    double p_case = pheno( 1 ) / ( pheno( 0 ) + pheno( 1 ) );

    *is_valid = arma::min( pheno ) >= 5;
    return log_double::from_log( pheno( 1 ) * log( p_case ) + pheno( 0 ) * log( 1 - p_case ) );
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
    
    arma::mat counts = joint_count( snp1, snp2, phenotype, weight ) + 0.5;
    arma::mat snp_pheno = zeros<mat>( 3, 2 );
    for(int i = 0; i < 3; i++)
    {
        snp_pheno( i, 0 ) = counts( 3*i, 0 ) + counts( 3*i + 1, 0 ) + counts( 3*i + 2, 0 );
        snp_pheno( i, 1 ) = counts( 3*i, 1 ) + counts( 3*i + 1, 1 ) + counts( 3*i + 2, 1 );
    }
    
    arma::vec p_snp = snp_pheno.col( 1 ) / ( snp_pheno.col( 0 ) + snp_pheno.col( 1 ) );
    
    *is_valid = arma::min( arma::min( snp_pheno ) ) >= 5;
    return log_double::from_log( accu( snp_pheno.col( 1 ) % arma::log( p_snp ) + snp_pheno.col( 0 ) % arma::log( 1 - p_snp ) ) );
}
