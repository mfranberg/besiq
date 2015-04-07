#include <bayesic/stats/loglinear_models.hpp>

using namespace arma;

full::full()
: loglinear_model::loglinear_model( 16, 0 )
{

}

log_double
full::prob(const arma::mat &count)
{
    arma::vec p_full = count.col( 1 ) / ( count.col( 0 ) + count.col( 1 ) );

    return log_double::from_log( accu( count.col( 1 ) % arma::log( p_full ) + count.col( 0 ) % arma::log( 1 - p_full ) ) );
}

block::block()
: loglinear_model::loglinear_model( 8, 8 )
{

}

log_double
block::prob(const arma::mat &count)
{
    arma::rowvec pheno = sum( count, 0 );
    double p_case = pheno( 1 ) / ( pheno( 0 ) + pheno( 1 ) );

    return log_double::from_log( pheno( 1 ) * log( p_case ) + pheno( 0 ) * log( 1 - p_case ) );
}

partial::partial(bool is_first)
: loglinear_model::loglinear_model( 10, 6 ),
  m_is_first( is_first )
{

}

log_double
partial::prob(const arma::mat &count)
{
    arma::mat snp_pheno = zeros<mat>( 3, 2 );
    for(int i = 0; i < 3; i++)
    {
        if( m_is_first )
        {
            snp_pheno( i, 0 ) = count( 3*i, 0 ) + count( 3*i + 1, 0 ) + count( 3*i + 2, 0 );
            snp_pheno( i, 1 ) = count( 3*i, 1 ) + count( 3*i + 1, 1 ) + count( 3*i + 2, 1 );
        }
        else
        {
            snp_pheno( i, 0 ) = count( 3*0 + i, 0 ) + count( 3*1 + i, 0 ) + count( 3*2 + i, 0 );
            snp_pheno( i, 1 ) = count( 3*0 + i, 1 ) + count( 3*1 + i, 1 ) + count( 3*2 + i, 1 );
        }
    }
    
    arma::vec p_snp = snp_pheno.col( 1 ) / ( snp_pheno.col( 0 ) + snp_pheno.col( 1 ) );
    
    return log_double::from_log( accu( snp_pheno.col( 1 ) % arma::log( p_snp ) + snp_pheno.col( 0 ) % arma::log( 1 - p_snp ) ) );
}
