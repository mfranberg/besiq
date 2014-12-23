#include <dcdflib/libdcdf.hpp>
#include <bayesic/method/wald_method.hpp>
#include <bayesic/stats/snp_count.hpp>

wald_method::wald_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_weight = arma::ones<arma::vec>( data->phenotype.n_elem );
}

void
wald_method::init(std::ostream &output)
{
    output << "LR\tP";
}

void
wald_method::run(const snp_row &row1, const snp_row &row2, std::ostream &output)
{
    arma::mat n = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    if( arma::min( arma::min( n ) ) <= 10 )
    {
        output << "NA\tNA";
        return;
    }

    arma::vec log_or0( 4 );
    arma::vec log_or1( 4 );
    arma::mat I0( 4, 4 );
    arma::mat I1( 4, 4 );
    int elem = 0;
    for(int i = 1; i <= 2; i++)
    {
        for(int j = 1; j <= 2; j++)
        {
            log_or0[ elem ] = log( n( i*3 + j, 0 ) * n( 0*3 + 0, 0 ) / n( 0*3 + j, 0 ) / n( i*3 + 0, 0 ) );
            log_or1[ elem ] = log( n( i*3 + j, 1 ) * n( 0*3 + 0, 1 ) / n( 0*3 + j, 1 ) / n( i*3 + 0, 1 ) );
            I0( elem, elem ) = 1.0 / n( i*3 + j, 0 ) + 1.0 / n( 0*3 + 0, 0 ) + 1.0 / n( 0*3 + j, 0 ) + 1.0 / n( i*3 + 0, 0 );
            I1( elem, elem ) = 1.0 / n( i*3 + j, 1 ) + 1.0 / n( 0*3 + 0, 1 ) + 1.0 / n( 0*3 + j, 1 ) + 1.0 / n( i*3 + 0, 1 );

            elem++;
        }
    }

    I0(0, 1) = 1.0 / n( 0*3 + 0, 0 ) + 1.0 / n( 1*3 + 0, 0 );
    I1(0, 1) = 1.0 / n( 0*3 + 0, 1 ) + 1.0 / n( 1*3 + 0, 1 );

    I0(0, 2) = 1.0 / n( 0*3 + 0, 0 ) + 1.0 / n( 0*3 + 1, 0 );
    I1(0, 2) = 1.0 / n( 0*3 + 0, 1 ) + 1.0 / n( 0*3 + 1, 1 ); 

    I0(0, 3) = 1.0 / n( 0*3 + 0, 0 );
    I1(0, 3) = 1.0 / n( 0*3 + 0, 1 ); 
    
    I0(1, 2) = I0(0, 3); 
    I1(1, 2) = I1(0, 3); 

    I0(1, 3) = 1.0 / n( 0*3 + 0, 0 ) + 1.0 / n( 0*3 + 2, 0 );
    I1(1, 3) = 1.0 / n( 0*3 + 0, 1 ) + 1.0 / n( 0*3 + 2, 1 ); 

    I0(2, 3) = 1.0 / n( 0*3 + 0, 0 ) + 1.0 / n( 2*3 + 0, 0 );
    I1(2, 3) = 1.0 / n( 0*3 + 0, 1 ) + 1.0 / n( 2*3 + 0, 1 );

    arma::vec log_diff = log_or0 - log_or1;
    arma::mat I = symmatu( I0 ) + symmatu( I1 );
    
    arma::mat Iinv( 4, 4 );
    if( !inv( Iinv, I ) )
    {
        output << "NA\tNA";
        return;
    }
    
    double chi = dot( log_diff, Iinv * log_diff );
    output << chi << "\t" << 1.0 - chi_square_cdf( chi, 4 );
}
