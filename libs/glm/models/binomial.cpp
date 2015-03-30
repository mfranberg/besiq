#include <glm/models/binomial.hpp>

#include <armadillo>

using namespace arma;

binomial::binomial(const std::string &link_name)
{
    m_link = make_link( link_name );
}

binomial::~binomial()
{
    delete m_link;
}

const glm_link &
binomial::get_link() const
{
    return *m_link;
}

bool
binomial::valid_mu(const arma::vec &mu) const
{
    for(int i = 0; i < mu.n_elem; i++)
    {
        if( mu[ i ] <= 0.0 || mu[ i ] >= 1.0 )
        {
            return false;
        }
    }
    
    return true;
}

vec
binomial::var(const arma::vec &mu) const
{
    return mu % ( 1.0 - mu );
}

double
binomial::likelihood(const arma::vec &mu, const arma::vec &y, const arma::uvec &missing) const
{
    double loglikelihood = 0.0;
    for(int i = 0; i < y.n_elem; i++)
    {
        if( missing[ i ] == 0 )
        {
            loglikelihood += y[ i ] * log( mu[ i ] ) + ( 1 - y[ i ] ) * log( 1 - mu[ i ] );
        }
    }

    return loglikelihood;
}
