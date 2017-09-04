#include "gene_environment.hpp"

gene_environment::gene_environment(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &phenotype, const std::vector<std::string> &cov_names)
        : m_genotypes( genotypes ),
          m_cov( cov ),
          m_phenotype( phenotype )
{
    
    std::vector<std::string> locus_names = genotypes->get_snp_names( );
    m_names.insert( m_names.end( ), locus_names.begin( ), locus_names.end( ) );
    if( cov_names.size( ) > 0 )
    {
        m_names.insert( m_names.end( ), cov_names.begin( ) + 2, cov_names.end( ) );
    }
}

void 
gene_environment::impute_missing()
{
    fill_missing_genotypes( );
    fill_missing_cov( );
    fill_missing_phenotypes( );
    compute_mean_sd( );
}

arma::vec
gene_environment::get_centered_phenotype()
{
    return m_phenotype - mean( m_phenotype );
}

size_t
gene_environment::get_num_samples()
{
    return m_phenotype.n_elem;
}

size_t
gene_environment::get_num_variables()
{
    return m_genotypes->get_snp_names( ).size( ) + m_cov.n_cols;
}

std::string
gene_environment::get_name(size_t index)
{
    return m_names[ index ];
}

std::vector<std::string> 
gene_environment::get_names()
{
    return m_names;
}

void
gene_environment::calculate_cor(const arma::vec &residual, arma::vec &c)
{
    int var = 0;
    #pragma omp parallel for
    for(int i = 0; i < m_genotypes->size( ); i++)
    {
        const snp_row &row = m_genotypes->get_row( i );
        double cor = 0.0;
        double m = m_mean[ i ];
        double s = m_sd[ i ];
        for(int j = 0; j < row.size( ); j++)
        {
            cor += ( ( row[ j ] - m ) / s ) * residual[ j ];
        }
        c[ i ] = cor;
    }

    var = m_genotypes->size( );
    for(int j = 0; j < m_cov.n_cols; j++)
    {
        double cor = dot( ( m_cov.col( j ) - m_mean[ var + j ] ) / m_sd[ var + j ],  residual );
        c[ var + j ] = cor;
    }
}

arma::vec
gene_environment::eig_prod(const arma::vec &u)
{
    arma::vec a = arma::zeros<arma::vec>( get_num_variables( ) );

    size_t var = 0;
    #pragma omp parallel for
    for(int i = 0; i < m_genotypes->size( ); i++)
    {
        const snp_row &row = m_genotypes->get_row( i );
        double m = m_mean[ i ];
        double s = m_sd[ i ];
        for(int j = 0; j < row.size( ); j++)
        {
            a[ i ] += ( ( row[ j ] - m ) / s ) * u[ j ];
        }
    }

    var = m_genotypes->size( );
    for(int i = 0; i < m_cov.n_cols; i++)
    {
        a[ var + i ] = dot( ( m_cov.col( i ) - m_mean[ var + i ] ) / m_sd[ var + i ], u );
    }

    return a;
}

arma::mat
gene_environment::get_active(const arma::uvec &active)
{
    size_t n = get_num_samples( );
    arma::mat X( n, active.n_elem );
    size_t var = 0;
    for(int i = 0; i < active.n_elem; i++)
    {
        double m = m_mean[ active[ i ] ];
        double s = m_sd[ active[ i ] ];
        if( active[ i ] < m_genotypes->size( ) )
        {
            const snp_row &row = m_genotypes->get_row( active[ i ] );
            for(int j = 0; j < row.size( ); j++)
            {
                X( j, var ) = ( row[ j ] - m ) / s;
            }
        }
        else
        {
            X.col( var ) = ( m_cov.col( active[ i ] - m_genotypes->size( ) ) - m ) / s;
        }
        var++;
    }

    return X;
}

arma::mat
gene_environment::get_active_raw(const arma::uvec &active)
{
    size_t n = get_num_samples( );
    arma::mat X( n, active.n_elem );
    size_t var = 0;
    for(int i = 0; i < active.n_elem; i++)
    {
        double m = m_mean[ active[ i ] ];
        double s = m_sd[ active[ i ] ];
        if( active[ i ] < m_genotypes->size( ) )
        {
            const snp_row &row = m_genotypes->get_row( active[ i ] );
            for(int j = 0; j < row.size( ); j++)
            {
                X( j, var ) = row[ j ];
            }
        }
        else
        {
            X.col( var ) = m_cov.col( active[ i ] - m_genotypes->size( ) );
        }
        var++;
    }

    return X;
}

void
gene_environment::fill_missing_genotypes()
{
    #pragma omp parallel for
    for(int i = 0; i < m_genotypes->size( ); i++)
    {
        snp_row &row = m_genotypes->get_row( i );
        int total = 0;
        double mu = 0;
        for(int j = 0; j < row.size( ); j++)
        {
            if( row[ j ] != 3 )
            {
                mu += row[ j ];
                total++;
            }
        }
        mu = mu / total;

        int mu_int = round( mu );
        for(int j = 0; j < row.size( ); j++)
        {
            if( row[ j ] == 3 )
            {
                row.assign( j, mu_int );
            }
        }
    }
}

void
gene_environment::fill_missing_phenotypes()
{
    double mu = 0.0;
    int total = 0;
    for(int i = 0; i < m_phenotype.n_elem; i++)
    {
        if( m_phenotype[ i ] == m_phenotype[ i ] )
        {
            mu += m_phenotype[ i ];
            total++;
        }
    }
    mu = mu / total;

    for(int i = 0; i < m_phenotype.n_elem; i++)
    {
        if( m_phenotype[ i ] != m_phenotype[ i ] )
        {
            m_phenotype[ i ] = mu;
        }
    }
}

void
gene_environment::fill_missing_cov()
{
    for(int i = 0; i < m_cov.n_cols; i++)
    {
        double mu = 0.0;
        int total = 0;
        for(int j = 0; j < m_cov.n_rows; j++)
        {
            // Check NaN, may fail in some compilers
            if( m_cov( j, i ) == m_cov( j, i ) )
            {
                mu += m_cov( j, i );
                total++;
            }
        }
        mu = mu / total;

        for(int j = 0; j < m_cov.n_rows; j++)
        {
            // Check NaN, may fail in some compilers
            if( m_cov( j, i ) != m_cov( j, i ) )
            {
                m_cov( j, i ) = mu;
            }
        }
    }
}

void
gene_environment::compute_mean_sd()
{
    m_mean = arma::zeros<arma::vec>( get_num_variables( ) );
    m_sd = arma::zeros<arma::vec>( get_num_variables( ) );
    size_t var = 0;
    size_t n = get_num_samples( );
    for(int i = 0; i < m_genotypes->size( ); i++)
    {
        snp_row &row = m_genotypes->get_row( i );
        double m = 0;
        int table[] = {0,0,0};
        for(int j = 0; j < row.size( ); j++)
        {
            m += row[ j ];
        }
        m = m / n;

        double s = 0;
        for(int j = 0; j < row.size( ); j++)
        {
            s += pow( row[ j ] - m, 2 );
        }

        m_mean[ var ] = m;
        m_sd[ var ] = sqrt( s );

        var++;
    }

    for(int i = 0; i < m_cov.n_cols; i++)
    {
        double m = arma::mean( m_cov.col( i ) );
        double s = arma::accu( pow( m_cov.col( i ) - m, 2 ) );

        m_mean[ var ] = m;
        m_sd[ var ] = sqrt( s );

        var++;
    }
}
