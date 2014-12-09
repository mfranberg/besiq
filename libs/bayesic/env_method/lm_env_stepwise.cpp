#include <bayesic/env_method/lm_env_stepwise.hpp>

#include <dcdflib/libdcdf.hpp>

lm_env_stepwise::lm_env_stepwise(method_data_ptr data, const arma::mat &E)
: method_env_type::method_env_type( data ),
  m_null_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 1 ),
  m_snp_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 1 + 2 ),
  m_env_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 1 + E.n_cols ),
  m_add_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 1 + 2 + E.n_cols ),
  m_alt_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 1 + 2 + E.n_cols + E.n_cols*2 ),
  m_E( E )
{
    /*
     * Null matrix.
     */
    m_null_matrix.col( 0 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    m_snp_matrix.col( 0 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    m_env_matrix.col( 0 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    m_add_matrix.col( 0 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    m_alt_matrix.col( 0 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    for(int i = 0; i < data->covariate_matrix.n_cols; i++)
    {
        m_null_matrix.col( i + 1 ) = data->covariate_matrix.col( i );
        m_snp_matrix.col( i + 1 + 2 ) = data->covariate_matrix.col( i );
        m_env_matrix.col( i + 1 + E.n_cols ) = data->covariate_matrix.col( i );
        m_add_matrix.col( i + 1 + 2 + E.n_cols ) = data->covariate_matrix.col( i );
        m_alt_matrix.col( i + 1 + 2 + E.n_cols + E.n_cols*2 ) = data->covariate_matrix.col( i );
    }

    m_env_matrix.cols( 1, 1 + m_E.n_cols - 1 ) = m_E;
    m_add_matrix.cols( 1 + 2, 1 + 2 + m_E.n_cols - 1 ) = m_E;
    m_alt_matrix.cols( 1 + 2, 1 + 2 + m_E.n_cols - 1 ) = m_E;
}

void
lm_env_stepwise::init(std::ostream &output)
{
    output << "P_null\tP_snp\tP_env\tP_add\t";
}

void
lm_env_stepwise::init_matrix_with_snp(const snp_row &row, arma::uvec &missing)
{
    for(int i = 0; i < row.size( ); i++)
    {
        if( row[ i ] == 3 )
        {
            missing[ i ] = 1;
            continue;
        }

        double s1 = row[ i ] == 1 ? 1.0 : 0.0;
        double s2 = row[ i ] == 2 ? 1.0 : 0.0;

        m_snp_matrix( i, 1 ) = s1;
        m_snp_matrix( i, 2 ) = s2;

        m_add_matrix( i, 1 ) = s1;
        m_add_matrix( i, 2 ) = s2;

        m_alt_matrix( i, 1 ) = s1;
        m_alt_matrix( i, 2 ) = s2;

        m_alt_matrix.row( i ).cols( 1 + 2 + m_E.n_cols, 1 + 2 + 2*m_E.n_cols - 1 ) = m_E.row( i ) * s1;
        m_alt_matrix.row( i ).cols( 1 + 2 + 2*m_E.n_cols, 1+ 2 + 3*m_E.n_cols - 1 ) = m_E.row( i ) * s2;
    }
}

void lm_env_stepwise::run(const snp_row &row, std::ostream &output)
{
    arma::uvec missing = get_data( )->missing;

    init_matrix_with_snp( row, missing );

    lm_info null_info;
    lm( m_null_matrix, get_data( )->phenotype, missing, null_info );
    
    lm_info snp_info;
    lm( m_snp_matrix, get_data( )->phenotype, missing, snp_info );
    
    lm_info env_info;
    lm( m_env_matrix, get_data( )->phenotype, missing, env_info );
    
    lm_info add_info;
    lm( m_add_matrix, get_data( )->phenotype, missing, add_info );

    lm_info alt_info;
    lm( m_alt_matrix, get_data( )->phenotype, missing, alt_info );

    if( null_info.success && snp_info.success && env_info.success && add_info.success && alt_info.success )
    {
        try
        {
            double LR_null = -2 *( null_info.logl - alt_info.logl );
            double p_null = 1.0 - chi_square_cdf( LR_null, m_alt_matrix.n_cols - m_null_matrix.n_cols );

            double LR_snp = -2 *( snp_info.logl - alt_info.logl );
            double p_snp = 1.0 - chi_square_cdf( LR_snp, m_alt_matrix.n_cols - m_snp_matrix.n_cols );

            double LR_env = -2 *( env_info.logl - alt_info.logl );
            double p_env = 1.0 - chi_square_cdf( LR_env, m_alt_matrix.n_cols - m_env_matrix.n_cols );

            double LR_add = -2 *( add_info.logl - alt_info.logl );
            double p_add = 1.0 - chi_square_cdf( LR_add, m_alt_matrix.n_cols - m_add_matrix.n_cols );

            output << p_null << "\t" << p_snp << "\t" << p_env << "\t" << p_add << "\t";
        }
        catch(bad_domain_value &e)
        {
            output << "NA\tNA\tNA\tNA\t";
        }
    }
    else
    {
        output << "NA\tNA\tNA\tNA\t";
    }
}
