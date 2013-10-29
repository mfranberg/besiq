#include <bayesic/method/bayesic_method.hpp>
#include <glm/models/binomial.hpp>
#include <glm/irls.hpp>

bayesic_method::bayesic_method(method_data_ptr data, arma::vec alpha)
: method_type::method_type( data )
{
    double int_prior = 1.0 / ( 4.0 * data->num_interactions );
    double single_prior = 1.0 / ( 2.0 * data->num_single );

    if( data->single_prior > 0.0 )
    {
        double p = data->single_prior;
        int_prior = p*p;
        single_prior = p*(1-p);
    }

    m_models.push_back( new saturated( int_prior, alpha ) );
    m_models.push_back( new ld_assoc( single_prior, alpha, true ) );
    m_models.push_back( new ld_assoc( single_prior, alpha, false ) );
    m_models.push_back( new null( 1.0 - int_prior - 2*single_prior, alpha ) );
}

bayesic_method::~bayesic_method()
{
    for(int i = 0; i < m_models.size( ); i++)
    {
        delete m_models[ i ];
    }

    m_models.clear( );
}

void
bayesic_method::set_models(const std::vector<model *> &models)
{
    for(int i = 0; i < m_models.size( ); i++)
    {
        delete m_models[ i ];
    }

    m_models.clear( );
    m_models.assign( models.begin( ), models.end( ) );
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
        arma::mat design_matrix = get_data( )->covariate_matrix;
        design_matrix.insert_cols( 0, arma::ones<arma::vec>( get_data( )->phenotype.n_elem ) );

        binomial model;
        irls_info full_info;
        irls( design_matrix, get_data( )->phenotype, get_data( )->missing, model, full_info );

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
