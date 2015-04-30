#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/stagewise_method.hpp>
#include <bayesic/stats/snp_count.hpp>
#include <bayesic/stats/binomial_models.hpp>
#include <bayesic/stats/normal_models.hpp>

stagewise_method::stagewise_method(method_data_ptr data, const std::string &model)
: method_type::method_type( data ),
  m_model( model )
{
    if( model == "normal" )
    {
        m_models.push_back( new normal_full( ) );
        m_models.push_back( new normal_null( ) );
        m_models.push_back( new normal_single( true ) );
        m_models.push_back( new normal_single( false ) );
    }
    else if( model == "binomial" )
    {
        m_models.push_back( new binomial_full( ) );
        m_models.push_back( new binomial_null( ) );
        m_models.push_back( new binomial_single( true ) );
        m_models.push_back( new binomial_single( false ) );
    }

    m_weight = 1.0 - arma::conv_to<arma::vec>::from( data->missing );
}

std::vector<std::string>
stagewise_method::init()
{
    std::vector<std::string> header;

    header.push_back( "P_null" );
    header.push_back( "P_snp1" );
    header.push_back( "P_snp2" );

    return header;
}

void
stagewise_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    std::vector<log_double> likelihood( m_models.size( ), 0.0 );

    arma::mat count;
    float min_samples = 0.0;
    unsigned int sample_threshold = METHOD_SMALLEST_CELL_SIZE_BINOMIAL;
    if( m_model == "binomial" )
    {
        count = joint_count( row1, row2, get_data( )->phenotype, m_weight );
        set_num_ok_samples( (size_t) arma::accu( count ) );
        min_samples = arma::min( arma::min( count ) );
    }
    else if( m_model == "normal" )
    {
        count = joint_count_cont( row1, row2, get_data( )->phenotype, m_weight );
        set_num_ok_samples( (size_t) arma::accu( count.col( 1 ) ) );
        min_samples = arma::min( count.col( 1 ) );
        sample_threshold = METHOD_SMALLEST_CELL_SIZE_NORMAL;
    }
    
    if( min_samples < sample_threshold )
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
            output[ i - 1 ] = 1.0 - chi_square_cdf( LR, m_models[ 0 ]->df( ) - m_models[ i ]->df( ) );
        }
        catch(bad_domain_value &e)
        {
        }
    }
}
