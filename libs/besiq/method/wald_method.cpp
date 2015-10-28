#include <dcdflib/libdcdf.hpp>
#include <besiq/method/wald_method.hpp>
#include <besiq/stats/snp_count.hpp>

wald_method::wald_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_weight = arma::ones<arma::vec>( data->phenotype.n_elem );
}

std::vector<std::string>
wald_method::init()
{
    std::vector<std::string> header;

    header.push_back( "LR" );
    header.push_back( "P" );
    header.push_back( "df" );

    return header;
}

arma::mat
wald_method::get_last_C()
{
    return m_C;
}

arma::vec
wald_method::get_last_beta()
{
    return m_beta;
}

double
wald_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat n0 = arma::zeros<arma::mat>( 3, 3 );
    arma::mat n1 = arma::zeros<arma::mat>( 3, 3 );
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] == 3 || row2[ i ] == 3 || get_data( )->missing[ i ] == 1 )
        {
            continue;
        }

        unsigned int pheno = get_data( )->phenotype[ i ];
        if( pheno == 0 )
        {
            n0( row1[ i ], row2[ i ] ) += 1;
        }
        else if( pheno == 1 )
        {
            n1( row1[ i ], row2[ i ] ) += 1;
        }
    }

    arma::mat eta( 3, 3 );
    double num_samples = 0.0;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if( n0( i, j ) < METHOD_SMALLEST_CELL_SIZE_BINOMIAL || n1( i, j ) < METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
            {
                continue;
            }

            eta( i, j ) = log( n1( i, j ) / n0( i, j ) );
            num_samples += n1( i, j ) + n0( i, j );
        }
    }

    /* Find valid parameters and estimate beta */
    int num_valid = 0;
    arma::uvec valid( 4 );
    m_beta = arma::zeros<arma::vec>( 4 );
    int i_map[] = { 1, 1, 2, 2 };
    int j_map[] = { 1, 2, 1, 2 };
    for(int i = 0; i < 4; i++)
    {
        int c_i = i_map[ i ];
        int c_j = j_map[ i ];
        if( n0( 0, 0 ) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n0( 0, c_j ) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n0( c_i, 0 ) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n0( c_i, c_j) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n1( 0, 0 ) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n1( 0, c_j ) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n1( c_i, 0 ) >= METHOD_SMALLEST_CELL_SIZE_NORMAL &&
            n1( c_i, c_j) >= METHOD_SMALLEST_CELL_SIZE_NORMAL )
        {
            valid[ num_valid ] = i;
            m_beta[ num_valid ] = eta( 0, 0 ) - eta( 0, c_j ) - eta( c_i, 0 ) + eta( c_i, c_j );
            num_valid++;
        }
    }
    set_num_ok_samples( (size_t)num_samples );
    if( num_valid <= 0 )
    {
        return -9;
    }
    
    valid.resize( num_valid );
    m_beta.resize( num_valid );

    /* Construct covariance matrix */
    m_C = arma::zeros<arma::mat>( num_valid, num_valid );
    for(int iv = 0; iv < num_valid; iv++)
    {
        int i = valid[ iv ];
        int c_i = i_map[ i ];
        int c_j = j_map[ i ];

        for(int jv = 0; jv < num_valid; jv++)
        {
            int j = valid[ jv ];
            int o_i = i_map[ j ];
            int o_j = j_map[ j ];

            int same_row = c_i == o_i;
            int same_col = c_j == o_j;
            int in_cell = i == j;

            m_C( iv, jv ) = 1.0 / n0( 0, 0 ) + same_col / n0( 0, c_j ) + same_row / n0( c_i, 0 ) + in_cell / n0( c_i, c_j );
            m_C( iv, jv ) += 1.0 / n1( 0, 0 ) + same_col / n1( 0, c_j ) + same_row / n1( c_i, 0 ) + in_cell / n1( c_i, c_j );
        }
    }

    arma::mat Cinv( num_valid, num_valid );
    if( !inv( Cinv, m_C ) )
    {
        return -9;
    }
    
    /* Test if b != 0 with Wald test */
    double chi = dot( m_beta, Cinv * m_beta );
    output[ 0 ] = chi;
    output[ 1 ] = 1.0 - chi_square_cdf( chi, num_valid );
    output[ 2 ] = valid.n_elem;

    return output[ 1 ];
}
