#include <dcdflib/libdcdf.hpp>
#include <glm/models/binomial.hpp>
#include <glm/models/normal.hpp>

#include <bayesic/method/scaleinv_method.hpp>

scaleinv_method::scaleinv_method(method_data_ptr data, model_matrix &model_matrix, bool is_lm)
: method_type::method_type( data ),
  m_model_matrix( model_matrix )
{
    if( !is_lm )
    {
        m_model.push_back( new binomial( "identity" ) );
        m_model.push_back( new binomial( "log" ) );
        m_model.push_back( new binomial( "logc" ) );
        m_model.push_back( new binomial( "odds" ) );
        m_model.push_back( new binomial( "logit" ) );
    }
    else
    {
        m_model.push_back( new normal( "identity" ) );
    }
}

scaleinv_method::~scaleinv_method()
{
    for(int i = 0; i < m_model.size( ); i++)
    {
        delete m_model[ i ];
    }
}

std::vector<std::string>
scaleinv_method::init()
{
    std::vector<std::string> header;
    for(int i = 0; i < m_model.size( ); i++)
    {
        header.push_back( m_model[ i ]->get_link( ).get_name( ) );
    }
    return header;
}

void scaleinv_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::uvec missing = get_data( )->missing;
    m_model_matrix.update_matrix( row1, row2, missing );
    set_num_ok_samples( missing.n_elem - sum( missing ) );
    
    for(int i = 0; i < m_model.size( ); i++)
    {
        glm_info alt_info;
        glm_fit( m_model_matrix.get_alt( ), get_data( )->phenotype, missing, *m_model[ i ], alt_info );

        glm_info null_info;
        glm_fit( m_model_matrix.get_null( ), get_data( )->phenotype, missing, *m_model[ i ], null_info );

        if( !null_info.success || !alt_info.success )
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
