#ifndef __DIRICHLET_H__
#define __DIRICHLET_H__

#include <armadillo>

#include <cmath>
#include <numeric>
#include <vector>
#include <tr1/random>

/**
 * Computes the dirichlet multinomial probability of a vector
 * x with prior parameter alpha.
 * 
 * @param x The observations.
 * @param alpha The prior parameters of the dirichlet density.
 *
 * @return The posterior probability of x.
 */
double dirmult(const arma::vec &x, const arma::vec &alpha);

/**
 * Computes the dirichlet multinomial log probability of a vector
 * x with prior parameter alpha.
 * 
 * @param x The observations.
 * @param alpha The prior parameters of the dirichlet density.
 *
 * @return The log posterior probability of x.
 */
double ldirmult(const arma::vec &x, const arma::vec &alpha);

/**
 * Computes the log of the binomial coefficient (n choose k).
 *
 * @param n Number of elements to draw from.
 * @param k Number of elements to draw.
 *
 * @return The log of the binomial coefficient.
 */
double lbinomial(double n, double k);

/**
 * This class is responsible for generating samples from a
 * dirichlet distribution. It does so by using the fact that
 * gamma(a_i, b) / sum_i gamma( a_i, b ) is dirichlet distributed
 * with parameter a_i.
 */
class dir_generator
{
public:
    /**
     * Constructor.
     *
     * Initializes the random generator with the given seed.
     *
     * @param seed The seed given to the random generator.
     */
    dir_generator(unsigned long seed);

    /**
     * Generates a random sample from the dirichlet distribution
     * with the given parameters.
     *
     * @param x The parameters of the dirichlet density.
     *
     * @return A sample from the dirichlet density.
     */
    arma::vec sample(const arma::vec &alpha);

private:
    /**
     * Mersenne twister random generator.
     */
    std::tr1::mt19937 m_generator;
};

#endif /* End of __DIRICHLET_H__ */
