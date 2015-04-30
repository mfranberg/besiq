#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/caseonly_method.hpp>
#include <bayesic/stats/snp_count.hpp>

caseonly_method::caseonly_method(method_data_ptr data, const std::string &method)
: method_type::method_type( data ),
  m_method( method ),
  m_wald( data )
{
    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

std::vector<std::string>
caseonly_method::init()
{
    std::vector<std::string> header;
    if( m_method == "r2" || m_method == "css" )
    {
        if( m_method == "r2" )
        {
            header.push_back( "R2" );
        }
        else if( m_method == "css" )
        {
            header.push_back( "CSS" );
        }
        
        header.push_back( "P_ld" );
        std::vector<std::string> wald_header = m_wald.init( );
        header.insert( header.end( ), wald_header.begin( ), wald_header.end( ) );
    }
    else if( m_method == "contrast" )
    {
        header.push_back( "LDdiff" );
        header.push_back( "P" );
    }

    return header;
}

void
caseonly_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    if( m_method == "r2" )
    {
        compute_r2( row1, row2, output );
        m_wald.run( row1, row2, &output[ 2 ] );
    }
    else if( m_method == "css" )
    {
        compute_css( row1, row2, output );
        m_wald.run( row1, row2, &output[ 2 ] );
    }
    else if( m_method == "contrast" )
    {
        compute_contrast( row1, row2, output );
    }
}

void
caseonly_method::compute_r2(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat counts = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    arma::vec snp_snp = sum( counts, 1 );
    arma::vec snp1 = arma::zeros<arma::vec>( 3 );
    arma::vec snp2 = arma::zeros<arma::vec>( 3 );
    arma::vec u = arma::zeros<arma::vec>( 3 );
    arma::vec v = arma::zeros<arma::vec>( 3 );
    double N = arma::accu( counts );
    set_num_ok_samples( (size_t) N );
    if( arma::min( arma::min( counts ) ) < METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        return;
    }

    for(int i = 0; i < 3; i++)
    {
        u[ i ] = i;
        v[ i ] = i;
        for(int j = 0; j < 3; j++)
        {
            snp1[ i ] += snp_snp[ 3 * i + j ];
            snp2[ j ] += snp_snp[ 3 * i + j ];
        }
    }
    
    double m = ( dot( u, snp1 ) * dot( v, snp2 ) ) / N;
    double R2 = 0.0;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            double term = u[ i ] * v[ j ] * snp_snp[ 3 * i + j ] - m;
            R2 += term * term;
        }
    }

    double s = ( dot( u % u, snp1 ) * dot( v % v, snp2 ) );
    R2 = R2 / s;
    
    output[ 0 ] = R2;
    output[ 1 ] = 1.0 - chi_square_cdf( R2, 1 );
}

void
caseonly_method::compute_css(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat counts = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    arma::vec snp_snp = sum( counts, 1 );
    arma::vec snp1 = arma::zeros<arma::vec>( 3 );
    arma::vec snp2 = arma::zeros<arma::vec>( 3 );
    if( arma::min( arma::min( counts ) ) < METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        return;
    }

    double N = arma::accu( counts );
    set_num_ok_samples( (size_t) N );

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            snp1[ i ] += snp_snp[ 3 * i + j ];
            snp2[ j ] += snp_snp[ 3 * i + j ];
        }
    }
    
    arma::vec f1 = snp1 / N;
    arma::vec f2 = snp2 / N;

    arma::vec o = arma::zeros<arma::vec>( 4 );
    o[ 0 ] = snp_snp[ 3 * 0 + 0 ] + snp_snp[ 3 * 1 + 0 ] + snp_snp[ 3 * 2 + 0 ] + snp_snp[ 3 * 0 + 1 ] + snp_snp[ 3 * 0 + 2 ];
    o[ 1 ] = snp_snp[ 3 * 1 + 1 ];
    o[ 2 ] = snp_snp[ 3 * 2 + 1 ] + snp_snp[ 3 * 1 + 2 ];
    o[ 3 ] = snp_snp[ 3 * 2 + 2 ];

    arma::vec e = arma::zeros<arma::vec>( 4 );
    e[ 0 ] = f1[ 0 ] * f2[ 0 ] + f1[ 1 ] * f2[ 0 ] + f1[ 2 ] * f2[ 0 ] + f1[ 0 ] * f2[ 1 ] + f1[ 0 ] * f2[ 2 ];
    e[ 1 ] = f1[ 1 ] * f2[ 1 ];
    e[ 2 ] = f1[ 1 ] * f2[ 2 ] + f1[ 2 ] * f2[ 1 ];
    e[ 3 ] = f1[ 2 ] * f2[ 2 ];

    double chi2 = sum( pow( o - N * e, 2 ) / ( N * e ) );
    output[ 0 ] = chi2;
    output[ 1 ] = 1.0 - chi_square_cdf( chi2, 3 );
}

void
caseonly_method::compute_contrast(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat counts = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    arma::vec snp_snp = sum( counts, 1 );

    if( arma::min( arma::min( counts ) ) < METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        return;
    }
    
    double N = arma::accu( counts );
    double N_controls = sum( counts.col( 0 ) );
    double N_cases = sum( counts.col( 1 ) );
    set_num_ok_samples( (size_t) N );

    double p_a = (2 * (snp_snp[ 0 ] + snp_snp[ 1 ] + snp_snp[ 2 ] ) + snp_snp[ 3 ] + snp_snp[ 4 ] + snp_snp[ 5 ] ) / ( 2 * N );
    double p_b = (2 * (snp_snp[ 0 ] + snp_snp[ 3 ] + snp_snp[ 6 ] ) + snp_snp[ 1 ] + snp_snp[ 4 ] + snp_snp[ 7 ] ) / ( 2 * N );

    double p_ab_case = ( 0.5 * counts( 4, 1 ) + counts( 5, 1 ) + counts( 7, 1 ) + 2 * counts( 8, 1 ) ) / N_cases;
    double p_ab_control = ( 0.5 * counts( 4, 0 ) + counts( 5, 0 ) + counts( 7, 0 ) + 2 * counts( 8, 0 ) ) / N_controls;

    double delta_case = ( p_ab_case - p_a * p_b );
    double delta_control = ( p_ab_control - p_a * p_b );

    double sigma_case = sqrt( (p_a *(1-p_a) * p_b *(1-p_b )) / N_cases ); 
    double sigma_control = sqrt( (p_a *(1-p_a) * p_b *(1-p_b )) / N_controls );

    double z = (delta_case - delta_control) / sqrt( pow( sigma_case, 2 ) + pow( sigma_control, 2 ) );
    double cdf = norm_cdf( z, 0.0, 1.0 );

    output[ 0 ] = z;
    output[ 1 ] = 2 * std::min( cdf, 1 - cdf );
}
