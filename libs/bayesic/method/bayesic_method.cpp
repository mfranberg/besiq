#include <method/bayesic_method.hpp>
#include <models/binomial.hpp>
#include <irls.hpp>

bayesic_method::bayesic_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_models.push_back( new saturated( 1.0 / ( 4.0 * data->num_interactions ) ) );
    m_models.push_back( new ld_assoc( true, 1.0 / ( 4.0 * data->num_interactions ) ) );
    m_models.push_back( new ld_assoc( false, 1.0 / ( 4.0 * data->num_interactions ) ) );
    m_models.push_back( new null( 1.0 - ( 3.0 / ( 4.0 * data->num_interactions ) ) ) );
}

void
bayesic_method::init(std::ostream &output)
{
    output << "Posterior";
   
    if( get_data( )->covariate_matrix.n_elem == 0 )
    {
        m_weight = arma::ones<arma::vec>( get_data( )->missing.n_elem );
    }
    else
    {
        binomial model;
        irls_info full_info;
        irls( get_data( )->covariate_matrix, get_data( )->phenotype, get_data( )->missing, model, full_info );

        arma::mat null = arma::ones<arma::mat>( get_data( )->phenotype.size( ), 1 );
        irls_info null_info;
        irls( null, get_data( )->phenotype, get_data( )->missing, model, null_info );

        m_weight = abs( get_data( )->phenotype - full_info.mu ) / abs( get_data( )->phenotype - null_info.mu );
    }

    const arma::uvec &missing = get_data( )->missing;
    for(int i = 0; i < missing.size( ); i++)
    {
        if( missing[ i ] == 1 )
        {
            m_weight[ i ] = 0.0;
        }
    }
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
