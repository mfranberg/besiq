#include <dcdflib/libdcdf.hpp>
#include <glm/models/binomial.hpp>
#include <glm/models/normal.hpp>

#include <besiq/method/separate_method.hpp>

separate_method::separate_method(method_data_ptr data, glm_model *model)
: method_type::method_type( data ),
  m_model( model )
{
    m_model_matrix.push_back( new separate_matrix( data->covariate_matrix, data->phenotype.n_elem, DOM_DOM ) );
    m_model_matrix.push_back( new separate_matrix( data->covariate_matrix, data->phenotype.n_elem, REC_DOM ) );
    m_model_matrix.push_back( new separate_matrix( data->covariate_matrix, data->phenotype.n_elem, DOM_REC ) );
    m_model_matrix.push_back( new separate_matrix( data->covariate_matrix, data->phenotype.n_elem, REC_REC ) );
}

separate_method::~separate_method()
{
    for(int i = 0; i < m_model_matrix.size( ); i++)
    {
        delete m_model_matrix[ i ];
    }
}

std::vector<std::string>
separate_method::init()
{
    std::vector<std::string> header;
    header.push_back( "b_dd" );
    header.push_back( "P_dd" );
    header.push_back( "b_rd" );
    header.push_back( "P_rd" );
    header.push_back( "b_dr" );
    header.push_back( "P_dr" );
    header.push_back( "b_rr" );
    header.push_back( "P_rr" );
    
    return header;
}

void separate_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    size_t num_samples = get_data( )->missing.n_elem - sum( get_data( )->missing );
    
    for(int i = 0; i < m_model_matrix.size( ); i++)
    {
        arma::uvec missing = get_data( )->missing;
        m_model_matrix[ i ]->update_matrix( row1, row2, missing );
        glm_info alt_info;
        arma::vec b = glm_fit( m_model_matrix[ i ]->get_alt( ), get_data( )->phenotype, missing, *m_model, alt_info );

        glm_info null_info;
        glm_fit( m_model_matrix[ i ]->get_null( ), get_data( )->phenotype, missing, *m_model, null_info );
        num_samples = std::min( (size_t) (missing.n_elem - sum( missing )), num_samples );

        if( !null_info.success || !alt_info.success )
        {
            continue;
        }

        try
        {
            double LR = -2 * ( null_info.logl - alt_info.logl );
            double p = 1.0 - chi_square_cdf( LR, m_model_matrix[ i ]->num_df( ) );
            output[ 2*i ] = b[ 2 ];
            output[ 2*i + 1 ] = p;
        }
        catch(bad_domain_value &e)
        {
        }
    }
    
    set_num_ok_samples( num_samples );
}
