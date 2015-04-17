#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/loglinear_method.hpp>
#include <bayesic/stats/snp_count.hpp>

loglinear_method::loglinear_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_models.push_back( new binomial_full( ) );
    m_models.push_back( new binomial_single( true ) );
    m_models.push_back( new binomial_single( false ) );
    m_models.push_back( new binomial_null( ) );

    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

std::vector<std::string>
loglinear_method::init()
{
    std::vector<std::string> header;

    header.push_back( "P" );

    return header;
}

void
loglinear_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat count = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    size_t num_samples = arma::accu( count );
    set_num_ok_samples( num_samples );
    if( arma::min( arma::min( count ) ) < 10 )
    {
        return;
    }

    std::vector<log_double> likelihood( m_models.size( ), 0.0 );
    std::vector<double> bic( m_models.size( ), 0.0 );
    for(int i = 0; i < m_models.size( ); i++)
    {
        likelihood[ i ] = m_models[ i ]->prob( count );
        bic[ i ] = -2.0 * likelihood[ i ].log_value( ) + m_models[ i ]->df( ) * log( num_samples );
    }

    unsigned int best_model = std::distance( bic.begin( ), std::min_element( bic.begin( ) + 1, bic.end( ) ) );
    double LR = -2.0*(likelihood[ best_model ].log_value( ) - likelihood[ 0 ].log_value( ));

    try
    {
        double p_value = 1.0 - chi_square_cdf( LR, m_models[ 0 ]->df( ) - m_models[ best_model ]->df( ) );
        output[ 0 ] = p_value;
    }
    catch(bad_domain_value &e)
    {
    }
}
