#include <besiq/stats/snp_count.hpp>

using namespace arma;

arma::mat
joint_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    arma::mat counts = zeros<mat>( 9, 2 );
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 )
        {
            unsigned int pheno = (unsigned int) phenotype[ i ];
            counts( 3 * row1[ i ] + row2[ i ], pheno ) += weight[ i ];
        }
    }

    return counts;
}

arma::mat
joint_count_cont(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    arma::mat counts = zeros<mat>( 9, 3 );
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 )
        {
            double pheno = phenotype[ i ];
            counts( 3 * row1[ i ] + row2[ i ], 0 ) += weight[ i ] * pheno;
            counts( 3 * row1[ i ] + row2[ i ], 1 ) += weight[ i ] * 1.0;
            counts( 3 * row1[ i ] + row2[ i ], 2 ) += weight[ i ] * pheno * pheno;
        }
    }

    return counts;
}

arma::vec
joint_count(const snp_row &row1, const snp_row &row2)
{
    arma::vec counts = zeros<vec>( 9 );
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 )
        {
            /* XXX: Should we have a weight here? */
            counts[ 3 * row1[ i ] + row2[ i ] ] += 1.0;
        }
    }

    return counts;
}

arma::vec
pheno_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    arma::vec counts = zeros<vec>( 2 );
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 )
        {
            unsigned int pheno = (unsigned int) phenotype[ i ];
            counts[ pheno ] += weight[ i ];
        }
    }

    return counts;
}

arma::mat
single_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight)
{
    arma::mat counts = zeros<mat>( 3, 2 );
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 )
        {
            unsigned int pheno = (unsigned int) phenotype[ i ];
            counts( row1[ i ], pheno ) += weight[ i ];
        }
    }

    return counts;
}

arma::vec
compute_maf(const snp_row &row)
{
    double p = compute_real_maf( row );
    
    arma::vec counts = zeros<vec>( 3 );
    counts[ 0 ] = (1 - p)*(1 - p);
    counts[ 1 ] = 2 * p * ( 1 - p );
    counts[ 2 ] = p * p;

    return counts;
}

float
compute_real_maf(const snp_row &row)
{
    arma::vec counts = zeros<vec>( 3 );
    for(int i = 0; i < row.size( ); i++)
    {
        if( row[ i ] != 3 )
        {
            counts[ row[ i ] ] += 1.0;
        }
    }

    unsigned int total_count = sum( counts );
    if( total_count != 0 )
    {
        return ( counts[ 1 ] + 2.0 * counts[ 2 ] ) / ( 2.0 * total_count );
    }
    else
    {
        return 0.0;
    }
}
