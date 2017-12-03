#include "gene_environment.hpp"

gene_environment::gene_environment(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &phenotype, const std::vector<std::string> &cov_names, bool only_main)
    : m_genotypes( genotypes ),
      m_cov( cov ),
      m_phenotype( phenotype ),
      m_only_main( only_main )
{
    std::vector<std::string> locus_names = genotypes->get_snp_names( );
    m_names.insert( m_names.end( ), locus_names.begin( ), locus_names.end( ) );
    if( cov_names.size( ) > 0 )
    {
        m_names.insert( m_names.end( ), cov_names.begin( ) + 2, cov_names.end( ) );
    }

    if(!only_main)
    {
        for(int i = 2; i < cov_names.size(); i++)
        {
            for(int j = 0; j < locus_names.size(); j++)
            {
                m_names.push_back( cov_names[ i ] + ":" + locus_names[ j ] );
            }
        }
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
gene_environment::get_centered_phenotype() const
{
    return m_phenotype - mean( m_phenotype );
}

size_t
gene_environment::get_num_samples() const
{
    return m_phenotype.n_elem;
}

size_t
gene_environment::get_num_variables() const
{
    if( !m_only_main )
    {
        return m_genotypes->size( ) + m_cov.n_cols + m_genotypes->size( ) * m_cov.n_cols;
    }
    else
    {
        return m_genotypes->size( ) + m_cov.n_cols;
    }
}

std::string
gene_environment::get_name(size_t index) const
{
    return m_names[ index ];
}

std::vector<std::string> 
gene_environment::get_names() const
{
    return m_names;
}

void
gene_environment::calculate_cor(const arma::vec &residual, arma::vec &c) const
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
        double cor = dot( ( m_cov.col( j ) - m_mean[ var + j ] ) / m_sd[ var + j ], residual );
        c[ var + j ] = cor;
    }

    if( m_only_main )
    {
        return;
    }

    var = m_genotypes->size( ) + m_cov.n_cols;
    size_t n = get_num_samples( );
    #pragma omp parallel for
    for(int i = 0; i < m_cov.n_cols; i++)
    {
        double cov_mu = m_mean[ m_genotypes->size( ) + i ];
        for(int j = 0; j < m_genotypes->size( ); j++)
        {
            snp_row &row = m_genotypes->get_row( j );
            double geno_mu = m_mean[ j ];

            double cor = 0.0;
            double m = m_mean[ var + i * m_genotypes->size( ) + j ];
            double s = m_sd[ var + i * m_genotypes->size( ) + j ];
            for(int k = 0; k < n; k++)
            {
                cor += (( (m_cov( k, i ) - cov_mu) * (row[ k ] - geno_mu) - m ) / s ) * residual[ k ];
            }

            c[ var + i * m_genotypes->size( ) + j ] = cor;
        }
    }
}

arma::vec
gene_environment::eig_prod(const arma::vec &u) const
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

    if( m_only_main )
    {
        return a;
    }
    
    size_t n = get_num_samples( );
    var = m_genotypes->size( ) + m_cov.n_cols;
    #pragma omp parallel for
    for(int i = 0; i < m_cov.n_cols; i++)
    {
        double cov_mu = m_mean[ m_genotypes->size( ) + i ];
        for(int j = 0; j < m_genotypes->size( ); j++)
        {
            snp_row &row = m_genotypes->get_row( j );
            double geno_mu = m_mean[ j ];

            double prod = 0.0;
            double m = m_mean[ var + i * m_genotypes->size( ) + j ];
            double s = m_sd[ var + i * m_genotypes->size( ) + j ];
            for(int k = 0; k < n; k++)
            {
                prod += (( (m_cov( k, i ) - cov_mu) * (row[ k ] - geno_mu) - m ) / s ) * u[ k ];
            }

            a[ var + i * m_genotypes->size( ) + j ] = prod;
        }
    }

    return a;
}

arma::mat
gene_environment::get_active(const arma::uvec &active) const
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
        else if( active[ i ] < m_genotypes->size( ) + m_cov.n_cols )
        {
            X.col( var ) = ( m_cov.col( active[ i ] - m_genotypes->size( ) ) - m ) / s;
        }
        else
        {
            unsigned int cov_index = (active[ i ] - m_genotypes->size( ) - m_cov.n_cols) / m_genotypes->size( );
            unsigned int snp_index = (active[ i ] - m_genotypes->size( ) - m_cov.n_cols) % m_genotypes->size( );
        
            double cov_mu = m_mean[ m_genotypes->size( ) + cov_index ];
            double geno_mu = m_mean[ snp_index ];

            const snp_row &row = m_genotypes->get_row( snp_index );
            const arma::vec &cov = m_cov.col( cov_index );
            for(int j = 0; j < n; j++)
            {
                X( j, var ) = ( (row[ j ] - geno_mu) * (cov[ j ] - cov_mu) - m ) / s;
            }
        }
        var++;
    }

    return X;
}

arma::mat
gene_environment::get_active_raw(const arma::uvec &active) const
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
        else if( active[ i ] < m_genotypes->size( ) + m_cov.n_cols )
        {
            X.col( var ) = m_cov.col( active[ i ] - m_genotypes->size( ) );
        }
        else
        {
            unsigned int cov_index = (active[ i ] - m_genotypes->size( ) - m_cov.n_cols) / m_genotypes->size( );
            unsigned int snp_index = (active[ i ] - m_genotypes->size( ) - m_cov.n_cols) % m_genotypes->size( );
        
            double cov_mu = m_mean[ m_genotypes->size( ) + cov_index ];
            double geno_mu = m_mean[ snp_index ];

            const snp_row &row = m_genotypes->get_row( snp_index );
            const arma::vec &cov = m_cov.col( cov_index );
            for(int j = 0; j < n; j++)
            {
                X( j, var ) = row[ j ] * cov[ j ];
            }
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
    size_t n = get_num_samples( );
    
    size_t var = 0;
    #pragma omp parallel for
    for(int i = 0; i < m_genotypes->size( ); i++)
    {
        snp_row &row = m_genotypes->get_row( i );
        double m = 0;
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

        m_mean[ i ] = m;
        m_sd[ i ] = sqrt( s );
    }

    var = m_genotypes->size( );
    for(int i = 0; i < m_cov.n_cols; i++)
    {
        double m = arma::mean( m_cov.col( i ) );
        double s = arma::accu( pow( m_cov.col( i ) - m, 2 ) );

        m_mean[ var + i ] = m;
        m_sd[ var + i ] = sqrt( s );
    }

    if( m_only_main )
    {
        return;
    }

    var = m_genotypes->size( ) + m_cov.n_cols;
    #pragma omp parallel for
    for(int i = 0; i < m_cov.n_cols; i++)
    {
        double cov_mu = m_mean[ m_genotypes->size( ) + i ];
        for(int j = 0; j < m_genotypes->size( ); j++)
        {
            snp_row &row = m_genotypes->get_row( j );
            double geno_mu = m_mean[ j ];

            double mean = 0.0;
            for(int k = 0; k < n; k++)
            {
                mean += (m_cov( k, i ) - cov_mu) * ( row[ k ] - geno_mu );
            }
            mean = mean / n;

            double s = 0.0;
            for(int k = 0; k < n; k++)
            {
                s += pow( (m_cov( k, i ) - cov_mu) * (row[ k ] - geno_mu) - mean, 2 );
            }

            m_mean[ var + m_genotypes->size( ) * i + j ] = mean;
            m_sd[ var + m_genotypes->size( ) * i + j ] = sqrt( s );
        }
    }
}
