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

unsigned int
stepwise_method::num_ok_samples(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype)
{
    return arma::accu( joint_count( row1, row2, get_data( )->phenotype, m_weight ) + 1.0 );
}

void
stepwise_method::init(std::ostream &output)
{
    output << "P_null\tP_snp1\tP_snp2\tLD_LR";
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
stepwise_method::run(const snp_row &row1, const snp_row &row2, std::ostream &output)
{
    std::vector<log_double> likelihood( m_models.size( ), 0.0 );
    std::vector<double> bic( m_models.size( ), 0.0 );
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
                output << 1.0 - chi_square_cdf( LR, m_models[ i ]->df( ) ) << "\t";
            }
            catch(bad_domain_value &e)
            {
                output << "NA\t";
            }
        }
    }
    else
    {
        for(int i = 1; i < m_models.size( ); i++)
        {
            output << "NA\t";
        }
    }
    
    output << compute_ld_lr( row1, row2 );
}
