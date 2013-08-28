#include <method/bayesic_method.hpp>

bayesic_method::bayesic_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_models.push_back( new saturated( 1.0 / ( 4.0 * data->num_interactions ) ) );
    m_models.push_back( new ld_assoc( true, 1.0 / ( 4.0 * data->num_interactions ) ) );
    m_models.push_back( new ld_assoc( false, 1.0 / ( 4.0 * data->num_interactions ) ) );
    m_models.push_back( new null( 1.0 - ( 3.0 / ( 4.0 * data->num_interactions ) ) ) );

    m_weight = arma::ones<arma::vec>( data->missing.n_elem );
    for(int i = 0; i < data->missing.n_elem; i++)
    {
        m_weight[ i ] = 1 - data->missing[ i ];
    }
}

void
bayesic_method::init(std::ostream &output)
{
    output << "Posterior";
}

void bayesic_method::run(const snp_row &row1, const snp_row &row2, std::ostream &output)
{
    log_double denominator = 0.0;
    std::vector<log_double> prior_likelihood( m_models.size( ), 0.0 );
    for(int i = 0; i < m_models.size( ); i++)
    {
        prior_likelihood[ i ] = m_models[ i ]->prior( ) * m_models[ i ]->prob( row1, row2, get_data( )->phenotype, m_weight );
        denominator += prior_likelihood[ i ];  
    }
    log_double posterior = prior_likelihood[ 0 ] / denominator;

    std::cout << posterior.value( );
}
