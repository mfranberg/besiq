#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/caseonly_method.hpp>
#include <bayesic/stats/snp_count.hpp>

caseonly_method::caseonly_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

unsigned int
caseonly_method::num_ok_samples(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype)
{
    return arma::accu( joint_count( row1, row2, get_data( )->phenotype, m_weight ) + 1.0 );
}

std::vector<std::string>
caseonly_method::init()
{
    std::vector<std::string> header;
    header.push_back( "R2" );
    header.push_back( "P" );

    return header;
}

void
caseonly_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat counts = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    arma::vec snp_snp = sum( counts, 1 );
    arma::vec snp1 = arma::zeros<arma::vec>( 3 );
    arma::vec snp2 = arma::zeros<arma::vec>( 3 );
    arma::vec u = arma::zeros<arma::vec>( 3 );
    arma::vec v = arma::zeros<arma::vec>( 3 );
    double N = arma::accu( counts );

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
