#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include <stats/log_scale.hpp>

TEST(log_scale_test, Init)
{
    log_double a = 4.6;
    ASSERT_NEAR( a.value( ), 4.6, 0.00001 );
}

TEST(log_scale_test, Add)
{
    log_double a = 0.3;
    log_double b = 0.5;
    log_double c = a + b;

    ASSERT_NEAR( c.value( ), 0.8, 0.00001 );

    c = a + 0.5;
    ASSERT_NEAR( c.value( ), 0.8, 0.00001 );
    
    c = 0.3 + b;
    ASSERT_NEAR( c.value( ), 0.8, 0.00001 );

    c += c;
    ASSERT_NEAR( c.value( ), 1.6, 0.00001 );

    c += 0.4;
    ASSERT_NEAR( c.value( ), 2.0, 0.00001 );


    a = log_double::from_log( -6542.0 );
    b = log_double::from_log( -5312.0 );
    c = a + b;
    ASSERT_NEAR( c.log_value( ), -5312.0, 0.00001 );
}

TEST(log_scale_test, Sub)
{
    log_double a = 0.5;
    log_double b = 0.3;
    log_double c = a - b;

    ASSERT_NEAR( c.value( ), 0.2, 0.00001 );

    c = a - 0.3;
    ASSERT_NEAR( c.value( ), 0.2, 0.00001 );
    
    c = 0.5 - b;
    ASSERT_NEAR( c.value( ), 0.2, 0.00001 );

    c -= c;
    ASSERT_NEAR( c.value( ), 0.0, 0.00001 );

    ASSERT_THROW( b - a, std::domain_error );
}

TEST(log_scale_test, Mul)
{
    log_double a = 0.5;
    log_double b = 0.3;
    log_double c = a * b;

    ASSERT_NEAR( c.value( ), 0.15, 0.00001 );

    c = a * 0.3;
    ASSERT_NEAR( c.value( ), 0.15, 0.00001 );
    
    c = 0.5 * b;
    ASSERT_NEAR( c.value( ), 0.15, 0.00001 );

    c = 0.5;
    c *= c;
    ASSERT_NEAR( c.value( ), 0.25, 0.00001 );
}

TEST(log_scale_test, Div)
{
    log_double a = 0.25;
    log_double b = 0.5;
    log_double c = a / b;

    ASSERT_NEAR( c.value( ), 0.5, 0.00001 );

    c = a / 0.5;
    ASSERT_NEAR( c.value( ), 0.5, 0.00001 );
    
    c = 0.25 / b;
    ASSERT_NEAR( c.value( ), 0.5, 0.00001 );

    c /= c;
    ASSERT_NEAR( c.value( ), 1.0, 0.00001 );
}

TEST(log_scale_test, Cmp)
{
    log_double a = 0.25;
    log_double b = 0.5;

    ASSERT_TRUE( a < b );
    ASSERT_FALSE( b < a );
    ASSERT_FALSE( a < a );
    
    ASSERT_TRUE( a <= b );
    ASSERT_FALSE( b <= a );
    ASSERT_FALSE( a >= b );
    ASSERT_TRUE( b >= a );

    a = b;
    ASSERT_TRUE( a <= b );
    ASSERT_TRUE( a >= b );
    ASSERT_TRUE( a == b );
    ASSERT_FALSE( a != b );

    a = 0.0;
    ASSERT_TRUE( a < b );
    ASSERT_FALSE( b < a );
    ASSERT_FALSE( a < a );    
}

TEST(log_scale_test, Zero)
{
    log_double a = 0.5;
    log_double c = a + 0.0;

    ASSERT_NEAR( c.value( ), 0.5, 0.00001 );
    c = 0.0 + a;
    ASSERT_NEAR( c.value( ), 0.5, 0.00001 );
    
    c = a - 0.0;
    ASSERT_NEAR( c.value( ), 0.5, 0.00001 );
    ASSERT_THROW( 0.0 - a, std::domain_error );
    
    c = a * 0.0;
    ASSERT_NEAR( c.value( ), 0.0, 0.00001 );
    c = 0.0 * a;
    ASSERT_NEAR( c.value( ), 0.0, 0.00001 );

    c = 0.0 / a ;
    ASSERT_NEAR( c.value( ), 0.0, 0.00001 );

    ASSERT_THROW( a / 0.0, std::domain_error );
}

TEST(log_scale_test, FromLog)
{
    log_double a = log_double::from_log( log( 2.5 ) );
    ASSERT_NEAR( a.value( ), 2.5, 0.00001 );
}

