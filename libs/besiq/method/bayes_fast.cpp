#include <besiq/method/bayes_fast.hpp>

bayes_fast_method::bayes_fast_method(method_data_ptr data, arma::vec alpha)
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

bayes_fast_method::~bayes_fast_method()
{
    for(int i = 0; i < m_models.size( ); i++)
    {
        delete m_models[ i ];
    }

    m_models.clear( );
}

void
bayes_fast_method::set_models(const std::vector<model *> &models)
{
    for(int i = 0; i < m_models.size( ); i++)
    {
        delete m_models[ i ];
    }

    m_models.clear( );
    m_models.assign( models.begin( ), models.end( ) );
}

std::vector<std::string>
bayes_fast_method::init()
{
    m_weight = arma::ones<arma::vec>( get_data( )->missing.n_elem );
   
    const arma::uvec &missing = get_data( )->missing;
    for(int i = 0; i < missing.size( ); i++)
    {
        if( missing[ i ] == 1 )
        {
            m_weight[ i ] = 0.0;
        }
    }

    std::vector<std::string> header;
    header.push_back( "Posterior" );
    return header;
}

double bayes_fast_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    log_double denominator = 0.0;
    std::vector<log_double> prior_likelihood( m_models.size( ), 0.0 );
    for(int i = 0; i < m_models.size( ); i++)
    {
        prior_likelihood[ i ] = m_models[ i ]->prior( ) * m_models[ i ]->prob( row1, row2, get_data( )->phenotype, m_weight );
        denominator += prior_likelihood[ i ];
    }
    log_double posterior = prior_likelihood[ 0 ] / denominator;

    output[ 0 ] = posterior.value( );
    
    return posterior.value( );
}
