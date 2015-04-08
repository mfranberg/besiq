#include <dcdflib/libdcdf.hpp>
#include <glm/models/binomial.hpp>
#include <glm/irls.hpp>

#include <bayesic/method/scaleinv_method.hpp>

scaleinv_method::scaleinv_method(method_data_ptr data, model_matrix &model_matrix, bool is_lm)
: method_type::method_type( data ),
  m_model_matrix( model_matrix ),
  m_is_lm( is_lm )
{
    if( !m_is_lm )
    {
        m_header.push_back( "identity" );
        m_header.push_back( "log" );
        m_header.push_back( "logc" );
        m_header.push_back( "odds" );
        m_header.push_back( "logit" );
    }
    else
    {
        m_header.push_back( "identity" );
    }
}

std::vector<std::string>
scaleinv_method::init()
{
    return m_header;
}

void scaleinv_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::uvec missing = get_data( )->missing;
    m_model_matrix.update_matrix( row1, row2, missing );
    set_num_ok_samples( missing.n_elem - sum( missing ) );
    
    if( !m_is_lm )
    {
        for(int i = 0; i < m_header.size( ); i++)
        {
            binomial model( m_header[ i ] );
            irls_info alt_info;
            irls( m_model_matrix.get_alt( ), get_data( )->phenotype, missing, model, alt_info );

            irls_info null_info;
            irls( m_model_matrix.get_null( ), get_data( )->phenotype, missing, model, null_info );

            if( !null_info.converged || !alt_info.converged )
            {
                continue;
            }

            try
            {
                double LR = -2 * ( null_info.logl - alt_info.logl );
                double p = 1.0 - chi_square_cdf( LR, m_model_matrix.num_df( ) );
                output[ i ] = p;
            }
            catch(bad_domain_value &e)
            {
            }
        }
    }
    else
    {
        lm_info null_info;
        lm( m_model_matrix.get_null( ), get_data( )->phenotype, missing, null_info );

        lm_info alt_info;
        arma::vec b = lm( m_model_matrix.get_alt( ), get_data( )->phenotype, missing, alt_info );
        if( !null_info.success || !alt_info.success )
        {
            return;
        }

        try
        {
            double LR = -2 * ( null_info.logl - alt_info.logl );
            double p = 1.0 - chi_square_cdf( LR, m_model_matrix.num_df( ) );
            output[ 0 ] = p;
        }
        catch(bad_domain_value &e)
        {
        }
    }
}
