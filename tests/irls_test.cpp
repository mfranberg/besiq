#include <armadillo>
#include <gtest/gtest.h>

#include <glm/irls.hpp>
#include <glm/models/binomial.hpp>

using namespace arma;

TEST(IRLSTest, Simple)
{
    double A_aux[] = { 1.0, 1.0, 1.0, 1.0,
                       -1.160, -0.655, 0.4156, -1.740 };

    double y_aux[] = { 0.1131, 0.4271, 0.9694, 0.0164 };

    mat A( A_aux, 4, 2 );
    vec y( y_aux, 4 );
    binomial binomial_model( "logit" );

    glm_info output;
    vec b = irls( A, y, binomial_model, output );

    ASSERT_NEAR( b[ 0 ], 2.0, 0.01 ); 
    ASSERT_NEAR( b[ 1 ], 3.5, 0.01 ); 
}
