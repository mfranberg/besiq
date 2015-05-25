#include <assert.h>

#include <besiq/prior.hpp>

#include <besiq/stats/snp_count.hpp>
#include <besiq/stats/beta.hpp>

/**
 * Permutes the given phenotype and returns a new
 * vector containing the permuted phenotype.
 *
 * @param phenotype A vector containing phenotypes.
 * @param weight Weight for each individual, will be changed
 *               after this call to ensure that individuals with
 *               weight 0 always are ignored.
 *
 * @return Permuted phenotypes.
 */
arma::vec
permute_phenotype(const arma::vec &phenotype, arma::vec &weight)
{
    arma::vec permuted_phenotype = phenotype;
    for(int i = phenotype.n_elem - 1; i > 0; --i)
    {
        int j = random( ) % ( i + 1 );

        std::swap( permuted_phenotype[ i ], permuted_phenotype[ j ] );
        std::swap( weight[ i ], weight[ j ] );
    }
    
    return permuted_phenotype;
}

/**
 * Estimates the risk in a random cell and returns the probability.
 *
 * @param snp1 First snp.
 * @param snp2 Second snp.
 * @param phenotype Phenotype.
 * @param weight Weight for each sample.
 *
 * @return The estimated risk in a random cell.
 */
double
sample_risk(const snp_row &snp1, const snp_row &snp2, const arma::vec &phenotype, arma::vec &weight)
{
    int cell = random( ) % 9;
    arma::vec permuted_phenotype = permute_phenotype( phenotype, weight );
    arma::mat counts = joint_count( snp1, snp2, permuted_phenotype, weight );

    return ( counts( cell, 1 ) + 1 ) / ( counts( cell, 0 ) + counts( cell, 1 ) + 2 );
}

arma::vec
estimate_prior_parameters(genotype_matrix_ptr genotypes, const arma::vec &phenotype, const arma::uvec &missing, int num_samples)
{
    assert( genotypes->size( ) > 1 );
    assert( num_samples > 1 );
    assert( sum( missing ) < phenotype.n_elem - 2.0 );

    arma::vec weight = arma::ones<arma::vec>( phenotype.n_elem );
    for(int i = 0; i < phenotype.n_elem; i++)
    {
        if( missing[ i ] == 1 )
        {
            weight[ i ] = 0;
        }
    }
    
    arma::vec samples = arma::zeros<arma::vec>( num_samples );
    for(int i = 0; i < num_samples; i++)
    {
        int snp1 = random( ) % genotypes->size( );
        int snp2 = snp1;
        while( snp1 == snp2 )
        {
            snp2 = random( ) % genotypes->size( );
        }

        const snp_row &snp1_row = genotypes->get_row( snp1 );
        const snp_row &snp2_row = genotypes->get_row( snp2 );

        samples[ i ] = sample_risk( snp1_row, snp2_row, phenotype, weight );
    }

    return mom_beta( samples );
}
