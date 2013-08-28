#include <libdcdf.hpp>
#include <method/loglinear_method.hpp>
#include <stats/snp_count.hpp>

loglinear_method::loglinear_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_models.push_back( new full( ) );
    m_models.push_back( new partial( true ) );
    m_models.push_back( new partial( false ) );
    m_models.push_back( new block( ) );

    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

unsigned int
loglinear_method::num_ok_samples(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype)
{
    return arma::accu( joint_count( row1, row2, get_data( )->phenotype, m_weight ) + 1.0 );
}

void
loglinear_method::init(std::ostream &output)
{
    output << "P";
}

void
loglinear_method::run(const snp_row &row1, const snp_row &row2, std::ostream &output)
{
    double num_samples = num_ok_samples( row1, row2, get_data( )->phenotype );
    std::vector<log_double> likelihood( m_models.size( ), 0.0 );
    std::vector<double> bic( m_models.size( ), 0.0 );
    for(int i = 0; i < m_models.size( ); i++)
    {
        likelihood[ i ] = m_models[ i ]->prob( row1, row2, get_data( )->phenotype, m_weight );
        bic[ i ] = -2.0 * likelihood[ i ].log_value( ) + m_models[ i ]->num_params( ) * log( num_samples );
    }

    unsigned int best_model = std::distance( bic.begin( ), std::min_element( bic.begin( ) + 1, bic.end( ) ) );
    double LR = -2.0*(likelihood[ best_model ].log_value( ) - likelihood[ 0 ].log_value( ));
    
    output << 1.0 - chi_square_cdf( LR, m_models[ best_model ]->df( ) );
}
