#include <dcdflib/libdcdf.hpp>
#include <besiq/method/wald_separate_method.hpp>
#include <besiq/stats/snp_count.hpp>

wald_separate_method::wald_separate_method(method_data_ptr data, bool is_lm)
: method_type::method_type( data )
{
    m_weight = arma::ones<arma::vec>( data->phenotype.n_elem );
    m_is_lm = is_lm;
}

std::vector<std::string>
wald_separate_method::init()
{
    std::vector<std::string> header;

    header.push_back( "b_dd" );
    header.push_back( "P_dd" );
    header.push_back( "b_rd" );
    header.push_back( "P_rd" );
    header.push_back( "b_dr" );
    header.push_back( "P_dr" );
    header.push_back( "b_rr" );
    header.push_back( "P_rr" );

    return header;
}

void
wald_separate_method::compute_lm(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat n = joint_count_cont( row1, row2, get_data( )->phenotype, m_weight );

    size_t num_samples = arma::accu( n.col( 1 ) );
    set_num_ok_samples( num_samples );
    
    /* Calculate residual and estimate sigma^2 */
    double residual_sum = 0.0;
    arma::mat mu = arma::zeros<arma::mat>( 3, 3 );
    for(int i = 0; i < 9; i++)
    {
        double deviance = ( n( i, 2 ) - n( i, 0 ) * n( i, 0 ) / n( i, 1 ) );
        residual_sum += deviance;
    }
    double sigma2 = residual_sum / ( num_samples - 9 );

    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 1, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 3, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 4, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL )
    {
        double b_dd = n( 0, 0 ) / n( 0, 1 ) - n( 1, 0 ) / n( 1, 1 ) - n( 3, 0 ) / n( 3, 1 ) + n( 4, 0 ) / n( 4, 1 );
        double var_dd = sigma2 * ( 1.0 / n( 0, 1 ) + 1.0 / n( 1, 1 ) + 1.0 / n( 3, 1 ) + 1.0 / n( 4, 1 ) );
        double w_dd = b_dd * b_dd / var_dd;
        output[ 0 ] = b_dd;
        output[ 1 ] = 1.0 - chi_square_cdf( w_dd, 1 );
    }
    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 2, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 3, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 5, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL )
    {
        double b_rd = n( 0, 0 ) / n( 0, 1 ) - n( 2, 0 ) / n( 2, 1 ) - n( 3, 0 ) / n( 3, 1 ) + n( 5, 0 ) / n( 5, 1 );
        double var_rd = sigma2 * ( 1.0 / n( 0, 1 ) + 1.0 / n( 2, 1 ) + 1.0 / n( 3, 1 ) + 1.0 / n( 5, 1 ) );
        double w_rd = b_rd * b_rd / var_rd;
        output[ 2 ] = b_rd;
        output[ 3 ] = 1.0 - chi_square_cdf( w_rd, 1 );
    }
    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 1, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 6, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 7, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL )
    {
        double b_dr = n( 0, 0 ) / n( 0, 1 ) - n( 1, 0 ) / n( 1, 1 ) - n( 6, 0 ) / n( 6, 1 ) + n( 7, 0 ) / n( 7, 1 );
        double var_dr = sigma2 * ( 1.0 / n( 0, 1 ) + 1.0 / n( 1, 1 ) + 1.0 / n( 6, 1 ) + 1.0 / n( 7, 1 ) );
        double w_dr = b_dr * b_dr / var_dr;
        output[ 4 ] = b_dr;
        output[ 5 ] = 1.0 - chi_square_cdf( w_dr, 1 );
    }
    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 2, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 6, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL &&
        n( 8, 1 ) > METHOD_SMALLEST_CELL_SIZE_NORMAL )
    {
        double b_rr = n( 0, 0 ) / n( 0, 1 ) - n( 2, 0 ) / n( 2, 1 ) - n( 6, 0 ) / n( 6, 1 ) + n( 8, 0 ) / n( 8, 1 );
        double var_rr = sigma2 * ( 1.0 / n( 0, 1 ) + 1.0 / n( 2, 1 ) + 1.0 / n( 6, 1 ) + 1.0 / n( 8, 1 ) );
        double w_rr = b_rr * b_rr / var_rr;
        output[ 6 ] = b_rr;
        output[ 7 ] = 1.0 - chi_square_cdf( w_rr, 1 );
    }
}

void
wald_separate_method::compute_binomial(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat n = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    set_num_ok_samples( (size_t) arma::accu( n ) );

    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 1, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 3, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 4, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 0, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 1, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 3, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 4, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        double b_dd = log( n( 0, 1 ) / n( 0, 0 ) ) - log( n( 1, 1 ) / n( 1, 0 ) ) - log( n( 3, 1 ) / n( 3, 0 ) ) + log( n( 4, 1 ) / n( 4, 0 ) );
        double var_dd = 1.0 / n( 0, 0 ) + 1.0 / n( 0, 1 ) + 1.0 / n( 1, 0 ) + 1.0 / n( 1, 1 ) + 1.0 / n( 3, 0 ) + 1.0 / n( 3, 1 ) + 1.0 / n( 4, 0 ) + 1.0 / n( 4, 1 );
        double w_dd = b_dd * b_dd / var_dd;
        output[ 0 ] = b_dd;
        output[ 1 ] = 1.0 - chi_square_cdf( w_dd, 1 );
    }
    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 2, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 3, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 5, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 0, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 2, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 3, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 5, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        double b_rd = log( n( 0, 1 ) / n( 0, 0 ) ) - log( n( 2, 1 ) / n( 2, 0 ) ) - log( n( 3, 1 ) / n( 3, 0 ) ) + log( n( 5, 1 ) / n( 5, 0 ) );
        double var_rd = 1.0 / n( 0, 0 ) + 1.0 / n( 0, 1 ) + 1.0 / n( 2, 0 ) + 1.0 / n( 2, 1 ) + 1.0 / n( 3, 0 ) + 1.0 / n( 3, 1 ) + 1.0 / n( 5, 0 ) + 1.0 / n( 5, 1 );
        double w_rd = b_rd * b_rd / var_rd;
        output[ 2 ] = b_rd;
        output[ 3 ] = 1.0 - chi_square_cdf( w_rd, 1 );
    }
    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 1, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 6, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 7, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 0, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 1, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 6, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 7, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        double b_dr = log( n( 0, 1 ) / n( 0, 0 ) ) - log( n( 1, 1 ) / n( 1, 0 ) ) - log( n( 6, 1 ) / n( 6, 0 ) ) + log( n( 7, 1 ) / n( 7, 0 ) );
        double var_dr = 1.0 / n( 0, 0 ) + 1.0 / n( 0, 1 ) + 1.0 / n( 1, 0 ) + 1.0 / n( 1, 1 ) + 1.0 / n( 6, 0 ) + 1.0 / n( 6, 1 ) + 1.0 / n( 7, 0 ) + 1.0 / n( 7, 1 );
        double w_dr = b_dr * b_dr / var_dr;
        output[ 4 ] = b_dr;
        output[ 5 ] = 1.0 - chi_square_cdf( w_dr, 1 );
    }
    if( n( 0, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 2, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 6, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 8, 1 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 0, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 2, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 6, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL &&
        n( 8, 0 ) > METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        double b_rr = log( n( 0, 1 ) / n( 0, 0 ) ) - log( n( 2, 1 ) / n( 2, 0 ) ) - log( n( 6, 1 ) / n( 6, 0 ) ) + log( n( 8, 1 ) / n( 8, 0 ) );
        double var_rr = 1.0 / n( 0, 0 ) + 1.0 / n( 0, 1 ) + 1.0 / n( 2, 0 ) + 1.0 / n( 2, 1 ) + 1.0 / n( 6, 0 ) + 1.0 / n( 6, 1 ) + 1.0 / n( 8, 0 ) + 1.0 / n( 8, 1 );
        double w_rr = b_rr * b_rr / var_rr;
        output[ 6 ] = b_rr;
        output[ 7 ] = 1.0 - chi_square_cdf( w_rr, 1 );
    }
}

double
wald_separate_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    if( m_is_lm )
    {
        compute_lm( row1, row2, output );
    }
    else
    {
        compute_binomial( row1, row2, output );
    }

    return min_na( min_na( output[ 1 ], output[ 3 ] ),
                   min_na( output[ 5 ], output[ 7 ] ) );
}
