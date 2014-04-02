#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include <plink/snp_row.hpp>

TEST(snp_row_test, test_assign)
{
    snp_row row;
    row.resize( 50 );

    for(int i = 0; i < 50; i++)
    {
        row.assign( i, i % 4 );
    }

    for(int i = 0; i < 50; i++)
    {
        ASSERT_EQ( row[ i ], i % 4 );
    }
}
