#include <gtest/gtest.h>

#include <armadillo>
#include <cmath>
#include <stdexcept>

#include <bayesic/stats/snp_count.hpp>

class snp_count_test
: public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        /**
         * 1/1 0/1
         * 0/0 2/0
         */
        row1.resize( 5 );
        row1.assign( 0, 0 );
        row1.assign( 1, 1 );
        row1.assign( 2, 1 );
        row1.assign( 3, 0 );
        row1.assign( 4, 0 );

        row2.resize( 5 );
        row2.assign( 0, 0 );
        row2.assign( 1, 1 );
        row2.assign( 2, 1 );
        row2.assign( 3, 0 );
        row2.assign( 4, 1 );

        phenotype = arma::zeros<arma::vec>( 5 );
        phenotype[ 0 ] = 0;
        phenotype[ 1 ] = 1;
        phenotype[ 2 ] = 1;
        phenotype[ 3 ] = 1;
        phenotype[ 4 ] = 0;

        weight = arma::ones<arma::vec>( 5 );
    }

    snp_row row1;
    snp_row row2;
    arma::vec phenotype;
    arma::vec weight;

};


TEST_F(snp_count_test, joint_count1)
{
    arma::vec count = joint_count( row1, row2 );

    ASSERT_NEAR( count[ 0 ], 2.0, 0.00001 );
    ASSERT_NEAR( count[ 1 ], 1.0, 0.00001 );
    ASSERT_NEAR( count[ 2 ], 0.0, 0.00001 );
    ASSERT_NEAR( count[ 3 ], 0.0, 0.00001 );
    ASSERT_NEAR( count[ 4 ], 2.0, 0.00001 );
    ASSERT_NEAR( count[ 5 ], 0.0, 0.00001 );
    ASSERT_NEAR( count[ 6 ], 0.0, 0.00001 );
    ASSERT_NEAR( count[ 7 ], 0.0, 0.00001 );
    ASSERT_NEAR( count[ 8 ], 0.0, 0.00001 );
}

TEST_F(snp_count_test, joint_count2)
{
    arma::mat count = joint_count( row1, row2, phenotype, weight );

    ASSERT_NEAR( count( 0, 0 ), 1.0, 0.00001 );
    ASSERT_NEAR( count( 0, 1 ), 1.0, 0.00001 );
    ASSERT_NEAR( count( 1, 0 ), 1.0, 0.00001 );
    ASSERT_NEAR( count( 1, 1 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 2, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 2, 1 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 3, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 3, 1 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 4, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 4, 1 ), 2.0, 0.00001 );
    ASSERT_NEAR( count( 5, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 5, 1 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 6, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 6, 1 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 7, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 7, 1 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 8, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 8, 1 ), 0.0, 0.00001 );
}

TEST_F(snp_count_test, pheno_count)
{
    arma::vec count = pheno_count( row1, row2, phenotype, weight );

    ASSERT_NEAR( count[ 0 ], 2.0, 0.00001 );
    ASSERT_NEAR( count[ 1 ], 3.0, 0.00001 );
}


TEST_F(snp_count_test, single_count)
{
    arma::mat count = single_count( row1, row2, phenotype, weight );

    ASSERT_NEAR( count( 0, 0 ), 2.0, 0.00001 );
    ASSERT_NEAR( count( 0, 1 ), 1.0, 0.00001 );
    ASSERT_NEAR( count( 1, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 1, 1 ), 2.0, 0.00001 );
    ASSERT_NEAR( count( 2, 0 ), 0.0, 0.00001 );
    ASSERT_NEAR( count( 2, 1 ), 0.0, 0.00001 );
}

TEST(snp_count_maf_test, compute_maf)
{
    snp_row row;
    row.resize( 8 );

    row.assign( 0, 0 );
    row.assign( 1, 0 );
    row.assign( 2, 0 );
    row.assign( 3, 1 );
    row.assign( 4, 1 );
    row.assign( 5, 2 );
    row.assign( 6, 2 );
    row.assign( 7, 2 );

    
    arma::vec maf = compute_maf( row );
    ASSERT_NEAR( maf[ 0 ], 0.25, 0.00001 );
    ASSERT_NEAR( maf[ 1 ], 0.5, 0.00001 );
    ASSERT_NEAR( maf[ 2 ], 0.25, 0.00001 );
}
