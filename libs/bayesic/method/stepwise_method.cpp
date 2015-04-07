#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/stepwise_method.hpp>
#include <bayesic/stats/snp_count.hpp>

stepwise_method::stepwise_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_models.push_back( new full( ) );
    m_models.push_back( new block( ) );
    m_models.push_back( new partial( true ) );
    m_models.push_back( new partial( false ) );

    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

std::vector<std::string>
stepwise_method::init()
{
    std::vector<std::string> header;

    header.push_back( "P_null" );
    header.push_back( "P_snp1" );
    header.push_back( "P_snp2" );

    return header;
}

double
stepwise_method::compute_ld_lr(const snp_row &row1, const snp_row &row2)
{
    arma::mat n = joint_count( row1, row2 ) + 1.0;
    double total = arma::accu( n );
    double p_a = (2 * (n[ 0 ] + n[ 1 ] + n[ 2 ] ) + n[ 3 ] + n[ 4 ] + n[ 5 ] ) / ( 2 * total );
    double p_b = (2 * (n[ 0 ] + n[ 3 ] + n[ 6 ] ) + n[ 1 ] + n[ 4 ] + n[ 7 ] ) / ( 2 * total );
    double p_AA = ( n[ 0 ] + n[ 1 ] + n[ 2 ] ) / total;
    double p_BB = ( n[ 0 ] + n[ 3 ] + n[ 6 ] ) / total;

    double delta_ab = (1.0 / total) * ( 2 * n[ 0 ] + n[ 1 ] + n[ 3 ] + 0.5 * n[ 4 ] ) - 2 * p_a * p_b;

    double pi_a = p_a * ( 1 - p_a );
    double pi_b = p_b * ( 1 - p_b );

    double D_a = p_AA - p_a * p_a;
    double D_b = p_BB - p_b * p_b;

    double LR = total * delta_ab * delta_ab / ( ( pi_a + D_a ) * ( pi_b + D_b ) );

    return LR;
}

void
stepwise_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    std::vector<log_double> likelihood( m_models.size( ), 0.0 );
    arma::mat count = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    if( arma::min( arma::min( count ) ) < 10 )
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
