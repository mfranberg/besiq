#include <armadillo>
#include <sstream>

#include <gtest/gtest.h>

#include <covariates.hpp>

using namespace arma;

TEST(CovariateTest, Simple)
{
    std::stringstream cov_file;
    cov_file << "FID IID gender age\n";
    cov_file << "1 1 0 15\n";
    cov_file << "2 2 NA 35\n";
    cov_file << "3 3 1 20\n";
    cov_file << "5 5 15 100\n";

    std::vector<std::string> order;
    order.push_back( "1" );
    order.push_back( "2" );
    order.push_back( "3" );
    order.push_back( "4" );

    arma::uvec missing = arma::zeros<arma::uvec>( 4 );

    arma::mat cov = parse_covariate_matrix( cov_file, missing, order );

    ASSERT_EQ( cov.n_rows, 4 );
    ASSERT_EQ( cov.n_cols, 2 );

    ASSERT_NEAR( cov( 0, 0 ), 0.0, 0.0001 );
    ASSERT_NEAR( cov( 0, 1 ), 15.0, 0.0001 );
    ASSERT_NEAR( cov( 1, 0 ), 0.0, 0.0001 );
    ASSERT_NEAR( cov( 1, 1 ), 35.0, 0.0001 );
    ASSERT_NEAR( cov( 2, 0 ), 1.0, 0.0001 );
    ASSERT_NEAR( cov( 2, 1 ), 20.0, 0.0001 );
    ASSERT_NEAR( cov( 3, 0 ), 0.0, 0.0001 );
    ASSERT_NEAR( cov( 3, 1 ), 0.0, 0.0001 );

    ASSERT_EQ( missing[ 0 ], 0 );
    ASSERT_EQ( missing[ 1 ], 1 );
    ASSERT_EQ( missing[ 2 ], 0 );
    ASSERT_EQ( missing[ 3 ], 1 );
}

TEST(CovariateTest, HeaderFail)
{
    std::stringstream cov_file;
    cov_file << "FID random gender age\n";
    cov_file << "1 1 0 15\n";

    std::vector<std::string> order( 1, "1" );
    arma::uvec missing = arma::zeros<arma::uvec>( 1 );

    ASSERT_THROW( parse_covariate_matrix( cov_file, missing, order ), std::runtime_error );
}

TEST(CovariateTest, ColumnFail)
{
    std::stringstream cov_file;
    cov_file << "FID IID gender age\n";
    cov_file << "1 1 0\n";

    std::vector<std::string> order( 1, "1" );
    arma::uvec missing = arma::zeros<arma::uvec>( 1 );

    ASSERT_THROW( parse_covariate_matrix( cov_file, missing, order ), std::runtime_error );
}

TEST(CovariateTest, ParseFail)
{
    std::stringstream cov_file;
    cov_file << "FID IID gender age\n";
    cov_file << "1 1 0 hej\n";

    std::vector<std::string> order( 1, "1" );
    arma::uvec missing = arma::zeros<arma::uvec>( 1 );

    ASSERT_THROW( parse_covariate_matrix( cov_file, missing, order ), std::runtime_error );
}
