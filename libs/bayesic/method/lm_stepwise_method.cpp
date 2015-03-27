#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/lm_stepwise_method.hpp>
#include <bayesic/stats/snp_count.hpp>

lm_stepwise_method::lm_stepwise_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_models.push_back( new lm_full( ) );
    m_models.push_back( new intercept( ) );
    m_models.push_back( new single( true ) );
    m_models.push_back( new single( false ) );

    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

unsigned int
lm_stepwise_method::num_ok_samples(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype)
{
    return arma::accu( joint_count( row1, row2, get_data( )->phenotype, m_weight ) + 1.0 );
}

std::vector<std::string>
lm_stepwise_method::init()
{
    std::vector<std::string> header;
    header.push_back( "P_null" );
    header.push_back( "P_snp1" );
    header.push_back( "P_snp2" );

    return header;
}

void
lm_stepwise_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    std::vector<log_double> likelihood( m_models.size( ), 0.0 );
    bool all_valid = true;
    
    for(int i = 0; i < m_models.size( ); i++)
    {
        bool is_valid = false;
        likelihood[ i ] = m_models[ i ]->prob( row1, row2, get_data( )->phenotype, m_weight, &is_valid );
        all_valid = is_valid && all_valid;
    }

    if( all_valid )
    {

        for(int i = 1; i < m_models.size( ); i++)
        {
            double LR = -2.0*(likelihood[ i ].log_value( ) - likelihood[ 0 ].log_value( ));
            double p_value = 1.0 - chi_square_cdf( LR, m_models[ i ]->df( ) );

            try
            {
                output[ i - 1 ] = 1.0 - chi_square_cdf( LR, m_models[ i ]->df( ) );
            }
            catch(bad_domain_value &e)
            {
            }
        }
    }
}
