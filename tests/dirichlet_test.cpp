#include <gtest/gtest.h>

#include <stats/dirichlet.hpp>

using namespace arma;

TEST(DirichletTest, Simple)
{
    double x_data[] = { 1.0, 1.0 };
    double alpha_data[] = { 1.0, 1.0 };
    vec x( x_data, 2 );
    vec alpha( alpha_data, 2 ); 

    ASSERT_NEAR( dirmult( x, alpha ), 0.1666666, 0.0001 );

    double x_data2[] = { 0.0, 0.0 };
    vec x2( x_data2, 2 );
    ASSERT_NEAR( dirmult( x2, alpha ), 1.0, 0.0001 );
}

TEST(DirichletTest, Advanced)
{
    double x_data[] = { 100.0, 220.0 };
    double alpha_data[] = { 1.0, 1.0 };
    vec x( x_data, 2 );
    vec alpha( alpha_data, 2 ); 

    ASSERT_NEAR( ldirmult( x, alpha ), -201.48, 0.01 );
}
