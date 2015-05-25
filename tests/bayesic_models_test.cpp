#include <gtest/gtest.h>

#include <plink/snp_row.hpp>
#include <besiq/stats/besiq_models.hpp>

using namespace arma;

TEST(BayesicModelsTest, SIndependent)
{
    double prior = 1.0;
    arma::vec alpha = arma::ones<arma::vec>( 2 );
    sindependent model( prior, alpha, 50000 );
    
    snp_row snp1;
    snp1.resize( 4 );
    snp1.assign( 0, 0 );
    snp1.assign( 1, 0 );
    snp1.assign( 2, 0 );
    snp1.assign( 3, 0 );

    snp_row snp2;
    snp2.resize( 4 );
    snp2.assign( 0, 0 );
    snp2.assign( 1, 0 );
    snp2.assign( 2, 0 );
    snp2.assign( 3, 0 );
    
    arma::vec phenotype = arma::zeros<arma::vec>( 4 );
    phenotype[ 0 ] = 1;

    arma::vec weight = arma::ones<arma::vec>( 4 );

    log_double likelihood = model.prob( snp1, snp2, phenotype, weight );
    ASSERT_NEAR( likelihood.value( ), 0.09, 0.001 );
}
