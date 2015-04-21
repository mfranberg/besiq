#include <dcdflib/libdcdf.hpp>
#include <glm/models/binomial.hpp>
#include <glm/models/normal.hpp>
#include <glm/models/links/power.hpp>
#include <glm/models/links/power_odds.hpp>

#include <bayesic/method/boxcox_method.hpp>

boxcox_method::boxcox_method(method_data_ptr data, model_matrix &model_matrix, bool is_lm, float lambda_start, float lambda_end, float lambda_step)
: method_type::method_type( data ),
  m_model_matrix( model_matrix )
{
    float min_pheno = arma::min( data->phenotype % (1 - data->missing) );
    m_lambda2 = 0.0;
    if( min_pheno < 0 )
    {
        m_lambda2 = -min_pheno + 1;
    }

    for(float lambda = lambda_start; lambda <= lambda_end; lambda += lambda_step)
    {
        m_lambda.push_back( lambda );
        if( is_lm )
        {
            m_model.push_back( new normal( new power_link( lambda, m_lambda2 ) ) );
        }
        else
        {
            m_model.push_back( new binomial( new power_odds_link( lambda ) ) );
        }
    }
}

boxcox_method::~boxcox_method()
{
    for(int i = 0; i < m_model.size( ); i++)
    {
        delete m_model[ i ];
    }
}

std::vector<std::string>
boxcox_method::init()
{
    std::vector<std::string> header;
    header.push_back( "lambda" );
    header.push_back( "lambda2" );
    header.push_back( "LR" );
    header.push_back( "P" );
    
    return header;
}

void boxcox_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::uvec missing = get_data( )->missing;
    m_model_matrix.update_matrix( row1, row2, missing );
    set_num_ok_samples( missing.n_elem - sum( missing ) );
    

    double max_logl = -DBL_MAX;
    int best_index = -1;
    
    for(int i = 0; i < m_model.size( ); i++)
    {
        glm_info null_info;
        glm_fit( m_model_matrix.get_null( ), get_data( )->phenotype, missing, *m_model[ i ], null_info );
        
        if( !null_info.success )
        {
            continue;
        }

        if( null_info.logl > max_logl )
        {
            max_logl = null_info.logl;
            best_index = i;
        }
    }
    
    if( best_index == -1 )
    {
        return;
    }

    /* Fit alternative model and test against best null */
    glm_info alt_info;
    glm_fit( m_model_matrix.get_alt( ), get_data( )->phenotype, missing, *m_model[ best_index ], alt_info );

    if( alt_info.success )
    {

        try
        {
            double LR = -2 * ( max_logl - alt_info.logl );
            double p = 1.0 - chi_square_cdf( LR, m_model_matrix.num_df( ) );
            
            output[ 0 ] = m_lambda[ best_index ];
            output[ 1 ] = m_lambda2;
            output[ 2 ] = LR;
            output[ 3 ] = p;
        }
        catch(bad_domain_value &e)
        {
        }
    }
}
