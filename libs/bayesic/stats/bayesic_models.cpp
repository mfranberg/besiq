#include <stats/bayesic_models.hpp>

using namespace arma;

saturated::saturated(log_double prior)
: model::model( prior )
{

}


log_double
saturated::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    mat counts = joint_count( row1, row2, phenotype, weight );
    
    log_double likelihood = 1.0;
    for(int i = 0; i < counts.n_rows; i++)
    {
        vec row_count = counts.row( i ).t( );
        vec alpha = ones<vec>( 2 );
        likelihood *= log_double::from_log( ldirmult( row_count, alpha ) );
    }

    return likelihood;
}

null::null(log_double prior)
: model::model( prior )
{

}

log_double
null::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    return pheno_prob( row1, row2, phenotype, weight );
}

log_double
null::snp_prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    vec counts = joint_count( row1, row2 );
    vec alpha = compute_alpha( row1, row2 );

    return log_double::from_log( ldirmult( counts, alpha ) );
}

log_double
null::pheno_prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    vec counts = pheno_count( row1, row2, phenotype, weight );
    vec alpha = ones<vec>( 2 );

    return log_double::from_log( ldirmult( counts, alpha ) );
}

vec
null::compute_alpha(const snp_row &row1, const snp_row &row2)
{
    vec maf1 = compute_maf( row1 );
    vec maf2 = compute_maf( row2 );

    vec alpha = zeros<vec>( 9 );
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            alpha[ 3 * i + j ] = maf1[ i ] * maf2[ j ];
        }
    }

    return alpha;
}

ld_assoc::ld_assoc(log_double prior, bool is_first)
: model::model( prior ),
  m_is_first( is_first )
{

}

log_double
ld_assoc::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    const snp_row &snp1 = m_is_first ? row1 : row2;
    const snp_row &snp2 = m_is_first ? row2 : row1;

    mat counts = single_count( snp1, snp2, phenotype, weight );

    log_double likelihood = 1.0;
    for(int i = 0; i < counts.n_rows; i++)
    {
        vec row_count = counts.row( i ).t( );
        vec alpha = ones<vec>( 2 );
        likelihood *= log_double::from_log( ldirmult( row_count, alpha ) );
    }

    return likelihood;
}

sindependent::sindependent(log_double prior, int num_mc_iterations)
: model::model( prior ),
  m_rdir( 0 ),
  m_num_mc_iterations( num_mc_iterations )
{

}

log_double
sindependent::prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    mat counts = joint_count( row1, row2, phenotype, weight );

    vec alpha = ones<vec>( 2 );
    log_double likelihood = 0.0;
    for(int k = 0; k < m_num_mc_iterations; k++)
    {
        double p[ 3 ];
        double q[ 3 ];
        for(int i = 0; i < 3; i++)
        {
            p[ i ] = m_rdir.sample( alpha )[ 0 ];
            q[ i ] = m_rdir.sample( alpha )[ 0 ];
        }

        log_double factors = 1.0;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                double n_ij0 = counts( 3 * i + j, 0 );
                double n_ij1 = counts( 3 * i + j, 1 );

                log_double fac1 = log_double::from_log( n_ij0 * log( p[ i ] * q[ j ] ) );
                log_double fac2 = log_double::from_log( n_ij1 * log( 1.0 - p[ i ] * q[ j ] ) );
                factors *= fac1 * fac2;
            }
        }

        likelihood += factors / ( (double) m_num_mc_iterations );
    }

    return likelihood;
}
