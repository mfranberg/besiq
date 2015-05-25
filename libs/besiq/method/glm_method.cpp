#include <besiq/method/glm_method.hpp>

#include <dcdflib/libdcdf.hpp>

glm_method::glm_method(method_data_ptr data, const glm_model &model, model_matrix &model_matrix)
: method_type::method_type( data ),
  m_model( model ),
  m_model_matrix( model_matrix )
{
}

std::vector<std::string>
glm_method::init()
{
    std::vector<std::string> header;
    header.push_back( "LR" );
    header.push_back( "P" );

    return header;
}

void glm_method::run(const snp_row &row1, const snp_row &row2, float *output)
{ 
    arma::uvec missing = get_data( )->missing;

    m_model_matrix.update_matrix( row1, row2, missing );

    glm_info null_info;
    arma::vec b1 = glm_fit( m_model_matrix.get_null( ), get_data( )->phenotype, missing, m_model, null_info );

    glm_info alt_info;
    arma::vec b = glm_fit( m_model_matrix.get_alt( ), get_data( )->phenotype, missing, m_model, alt_info );

    set_num_ok_samples( missing.n_elem - sum( missing ) );

    if( null_info.success && alt_info.success )
    {
        double LR = -2 * ( null_info.logl - alt_info.logl );

        try
        {
            output[ 0 ] = LR;
            output[ 1 ] = 1.0 - chi_square_cdf( LR, m_model_matrix.num_df( ) );
        }
        catch(bad_domain_value &e)
        {

        }
    }
}
