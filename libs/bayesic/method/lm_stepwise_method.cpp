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

    m_weight = 1.0 - arma::conv_to<arma::vec>::from( data->missing );
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

    arma::mat count = joint_count_cont( row1, row2, get_data( )->phenotype, m_weight );
    set_num_ok_samples( (size_t) arma::accu( count.col( 1 ) ) );

    bool enough_samples = arma::min( count.col( 1 ) ) >= 10;
    if( !enough_samples )
    {
        return;
    }
    
    for(int i = 0; i < m_models.size( ); i++)
    {
        likelihood[ i ] = m_models[ i ]->prob( count );
    }

    for(int i = 1; i < m_models.size( ); i++)
    {
        double LR = -2.0*(likelihood[ i ].log_value( ) - likelihood[ 0 ].log_value( ));

        try
        {
            output[ i - 1 ] = 1.0 - chi_square_cdf( LR, m_models[ i ]->df( ) );
        }
        catch(bad_domain_value &e)
        {
        }
    }
}
