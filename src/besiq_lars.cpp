#include <iostream>
#include <sstream>
#include <cfloat>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>
#include <besiq/io/covariates.hpp>
#include <besiq/stats/snp_count.hpp>

#include <plink/plink_file.hpp>
#include <dcdflib/libdcdf.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-lars [OPTIONS] plink_file";
const std::string DESCRIPTION = "Uses the LARS algorithm to find important loci.";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

/**
 * This class is responsible for holding information about a single change
 * in one of the parameters in the model. P-values is computed when a parameter
 * first enters the model.
 */
struct lars_path
{
    /**
     * Constructor.
     *
     * @param pass The step in which the parameter enters the model.
     * @param variable Name of the variable.
     * @param beta Value of the parameter.
     * @param T Covariance test statistic.
     * @param pvalue Possible p-value of the parameter.
     * @param lambda The value of lambda at the knot.
     * @param beta_sum Current sum of betas.
     * @param explained_variance Current explained variance of the current model.
     */
    lars_path(int pass, std::string variable, double beta, double T, double pvalue, double lambda, double beta_sum, double explained_var)
        : m_pass( pass ),
          m_variable( variable ),
          m_beta( beta ),
          m_T( T ),
          m_pvalue( pvalue ),
          m_lambda( lambda ),
          m_beta_sum( beta_sum ),
          m_explained_var( explained_var )
    {
    }

    /**
     * The step in which the variable enters the model.
     */
    int m_pass;

    /**
     * Name of the variable.
     */
    std::string m_variable;

    /**
     * Value of the parameter.
     */
    double m_beta;
    
    /**
     * Test statistic
     */
    double m_T;

    /**
     * Possible p-value of the parameter (=1 when not the first time it enters the model).
     */
    double m_pvalue;

    /**
     * Lambda at the knot.
     */
    double m_lambda;

    /**
     * Sum of absolute values of parameters currently in the model.
     */
    double m_beta_sum;

    /**
     * Explained variance of current model.
     */
    double m_explained_var;
};

/**
 * This class is responsible for keeping track which variables
 * are currently in the active and inactive sets, excluding those
 * that are currently in the ignored set.
 */
class active_set
{
public:
    /**
     * Constructor.
     *
     * @param size Total number of variables ([0...size-1]).
     */
    active_set(size_t size)
    {
        m_active = std::set<unsigned int>( );
        m_ignored = std::set<unsigned int>( );
        m_size = size;
    }

    /**
     * Add a variable to the active set. It is assumed that this
     * variable is not ignored.
     *
     * @param x Index of variable to add.
     */
    void add(unsigned int x)
    {
        m_active.insert( x );
    }

    /**
     * Ignores a parameters, ignored parameters are neither
     * active or inactive.
     *
     * @param x Index of variable to ignore.
     */
    void ignore(unsigned int x)
    {
        m_active.erase( x );
        m_ignored.insert( x );
    }

    /**
     * Move a variable from active to inactive.
     *
     * @param Index of the variable to move.
     */
    void drop(unsigned int x)
    {
        m_active.erase( x );
    }

    /**
     * Returns a vector of indicies of the currently active variables
     * suitable for the arma package.
     *
     * @return A vector of indicies of active variables.
     */
    arma::uvec get_active()
    {
        arma::uvec active = arma::zeros<arma::uvec>( m_active.size( ) );
        std::set<unsigned int>::const_iterator it;
        unsigned int index = 0;
        for(it = m_active.begin( ); it != m_active.end( ); ++it)
        {
            active[ index ] = *it;
            index++;      
        }

        return active;
    }

    /**
     * Returns a vector of indicies of the currently inactive variables
     * suitable for the arma package.
     *
     * @return A vector of indicies of inactive variables.
     */
    arma::uvec get_inactive()
    {
        arma::uvec inactive = arma::zeros<arma::uvec>( m_size - m_active.size( ) - m_ignored.size( ) );
        unsigned int index = 0;
        for(int i = 0; i < m_size; i++)
        {
            if( m_active.count( i ) > 0 || m_ignored.count( i ) > 0 )
            {
                continue;
            }

            inactive[ index ] = i;
            index++;
        }

        return inactive;
    }

    /**
     * Returns the number of active variables.
     *
     * @return the number of active variables.
     */
    size_t size()
    {
        return m_active.size( );
    }

private:
    /**
     * Total number of variables.
     */
    size_t m_size;

    /**
     * Set of active variables.
     */
    std::set<unsigned int> m_active;

    /**
     * Set of ignored variables.
     */
    std::set<unsigned int> m_ignored;
};

/**
 * This class is responsible for both genetic and environmental data, along
 * with computing certain properties of this data that relates to the LARS
 * algorithm. Primarily, to avoid passing and operating on each data type
 * separately within the algorithm.
 */
class gene_environment
{
public:
    /**
     * Constructor.
     *
     * @param genotypes A matrix of genotypes.
     * @param cov A matrix of covariates.
     * @param phenotype A vector of phenotypes.
     * @param cov_names Names of the covariates.
     */
    gene_environment(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &phenotype, const std::vector<std::string> &cov_names)
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

    /**
     * Imputes the missing genotypes, covariates and phenotypes.
     */
    void impute_missing()
    {
        fill_missing_genotypes( );
        fill_missing_cov( );
        fill_missing_phenotypes( );
        compute_mean_sd( );
    }

    /**
     * Centers a phenotype and returns it.
     *
     * @return Returns a centered phenotype.
     */
    arma::vec get_centered_phenotype()
    {
        return m_phenotype - mean( m_phenotype );
    }

    /**
     * Returns the number of samples.
     *
     * @return the number of samples.
     */
    size_t get_num_samples()
    {
        return m_phenotype.n_elem;
    }

    /**
     * Returns the number of variables.
     *
     * @return the number of variables.
     */
    size_t get_num_variables()
    {
        return m_genotypes->get_snp_names( ).size( ) + m_cov.n_cols;
    }
   
    /**
     * Returns the name of the variable with the given index.
     *
     * @param index Index of the variable whoose name to return.
     *
     * @return The name of the variable with the given index.
     */ 
    std::string get_name(size_t index)
    {
        return m_names[ index ];
    }

    /**
     * Calculates the correlation between each variable and
     * the given vector of residuals. The results are stored
     * in the supplied vector c.
     *
     * @param residual A vector of residuals.
     * @param c A vector of at least the size of the number of variables,
     *          to store the computed correlations in.
     */
    void
    calculate_cor(const arma::vec &residual, arma::vec &c)
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
    
    /**
     * Computes the vector a from the LARS paper from the
     * vector u.
     *
     * @param u The u-vector from the LARS paper.
     *
     * @return The a-vector from the LARS paper.
     */
    arma::vec eig_prod(const arma::vec &u)
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

    /**
     * Returns a matrix of standardized variables according to the
     * given vector of indicides of active variables.
     *
     * @param active  Indicies of active variables.
     *
     * @return Matrix of standardized active variables.
     */
    arma::mat get_active(const arma::uvec &active)
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

private:
    /**
     * Assigns missing genotypes with genotype mean.
     */
    void fill_missing_genotypes()
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

    /**
     * Assigns missing phenotypes with phenotype mean.
     */
    void fill_missing_phenotypes()
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

    /**
     * Assigns missing covariates with covariate mean.
     */
    void fill_missing_cov()
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

    /**
     * Computes the mean and standard deviation of each variable
     * and stores it in the m_mean and m_sd vectors. This is performed
     * in this way because it is too expensive to store the genotypes
     * as floats.
     */
    void compute_mean_sd()
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

private:
    /**
     * Matrix of genotypes.
     */
    genotype_matrix_ptr m_genotypes;

    /**
     * Matrix of covariates.
     */
    arma::mat m_cov;

    /**
     * Vector of phenotypes.
     */
    arma::vec m_phenotype;

    /**
     * Vector of mean values for each variable.
     */
    arma::vec m_mean;

    /**
     * Vector of standard deviation for each variable.
     */
    arma::vec m_sd;

    /**
     * Vector variable names.
     */
    std::vector<std::string> m_names;
};

arma::vec
nice_division(arma::vec &a, arma::vec &b)
{
    b.elem( find( b < DBL_MIN ) ).fill( DBL_MIN );
    arma::vec x = a / b;
    x.elem( find( x <= 0 ) ).fill( DBL_MAX );
    
    return x;
}

/**
 * Solves the lasso problem for a given lambda. Used for computing the
 * null model during significance testing. The algorithm is currently
 * a pathwise coordinate descent.
 *
 * @param X Covariates.
 * @param y Outcome.
 * @param lambda Shrinkage factor.
 * @param start Starting values for beta.
 * @max_num_iter Maximum number of iterations before giving up.
 */
arma::vec optimize_lars_gd(const arma::mat &X, const arma::mat &y, double lambda, const arma::vec &start, size_t max_num_iter = 50)
{
    arma::vec beta = start;
    arma::vec r = y - X * beta;
    double prev_rss = 0;
    double cur_rss = sum( r % r );
    int num_iter = 0;

    while( std::abs( prev_rss - cur_rss ) / prev_rss > 1e-20 && num_iter++ < max_num_iter )
    {
        prev_rss = cur_rss;
        for(int j = 0; j < X.n_cols; j++)
        {
            double beta_star = arma::dot( X.col( j ), r ) + beta[ j ];
            double update = std::abs( beta_star ) - lambda;
            update = (update > 0) ? update : 0;

            double new_beta = (beta_star > 0) ? update : -update;
            double beta_increase = new_beta - beta[ j ];

            beta[ j ] = new_beta;

            if( std::abs( beta_increase ) > 0 )
            {
                r = r - X.col( j ) * beta_increase;
            }
        }
        cur_rss = sum( r % r );
    }

    return beta;
}

/**
 * Finds the maximum correlation and returns its index.
 *
 * @param v A vector of correlations.
 * @param inactive A subset of indicies to consider.
 * @param max_index The index of the maximum value in v will be stored here.
 *
 * @reutrn The maximum value of v restricted to the indicies in inactive.
 */
double
find_max(const arma::vec &v, const arma::uvec &inactive, unsigned int *max_index)
{
    double max_value = 0;
    *max_index = v.n_elem;
    for(int i = 0; i < inactive.n_elem; i++)
    {
        unsigned int cur_index = inactive[ i ];
        double cur_value = v[ cur_index ];

        // Notice the order of these two statements
        *max_index = ( cur_value > max_value ) ? cur_index : *max_index;
        max_value = ( cur_value > max_value ) ? cur_value : max_value;
    }

    return max_value;
}

std::vector<lars_path>
lars(gene_environment &variable_set, size_t max_vars = 15, bool lasso = true, bool only_p = false, double threshold = 1.0)
{
    variable_set.impute_missing( );
    arma::vec phenotype = variable_set.get_centered_phenotype( );
    size_t n = variable_set.get_num_samples( );
    size_t m = variable_set.get_num_variables( );

    arma::vec mu = arma::zeros<arma::vec>( n );
    arma::vec beta = arma::zeros<arma::vec>( m );
    arma::vec c = arma::zeros<arma::vec>( m );
    std::vector<lars_path> path;

    double pheno_mean = sum( phenotype ) / n;
    double pheno_var = sum( pow( phenotype - pheno_mean, 2 ) ) / (n - 1);
    double base_cor = sum( phenotype * pheno_mean );

    active_set active( m );
    arma::uvec inactive;
    arma::uvec drop;
    double eps = 1e-10;
    double last_p = threshold;

    int i = 0;
    while( active.size( ) < std::min( max_vars, m ) && last_p <= threshold )
    {
        variable_set.calculate_cor( phenotype - mu, c );
        arma::vec cabs = abs( c );
        unsigned int max_index;
        double C = find_max( cabs, active.get_inactive( ), &max_index );

        if( C < eps )
        {
            break;
        }
       
        arma::uvec prev_active = active.get_active( ); 
        active.add( max_index );
        inactive = active.get_inactive( );

        /* Create the active matrix */
        arma::mat X_active = variable_set.get_active( active.get_active( ) );

        /* Equation 2.4 */
        arma::vec s = sign( c.elem( active.get_active( ) ) );

        /* Equation 2.5 */
        arma::mat G = X_active.t( ) * X_active;
        arma::mat Ginv;
        if( !inv( Ginv, G ) )
        {
            active.ignore( max_index );
            continue;
        }
        double A = 1.0 / sqrt( arma::accu( s.t( ) * Ginv * s ) );

        /* Equation 2.6 */
        arma::vec w = A * Ginv * s;
        arma::vec u = X_active * w;

        double gamma = C / A;
        if( active.size( ) < m )
        { 
            /* Equation 2.11 */
            arma::vec a = variable_set.eig_prod( u );

            /* Equation 2.13 */
            arma::vec cc = c.elem( inactive );
            arma::vec ac = a.elem( inactive );
            
            arma::vec a1 = C - cc;
            arma::vec b1 = A - ac;
            arma::vec gn = nice_division( a1, b1 );

            arma::vec a2 = C + cc;
            arma::vec b2 = A + ac;
            arma::vec gp = nice_division( a2, b2 );

            gamma = std::min( gn.min( ), gp.min( ) );
        }

        if( lasso )
        {
            drop.clear( );
            /* Equation 3.4 */
            arma::vec gammaj = -beta.elem( active.get_active( ) ) / w;

            /* Equation 3.5 */
            arma::uvec valid_elems = arma::find( gammaj > eps );
            if( valid_elems.n_elem > 0 )
            {
                double gammatilde = std::min( arma::min( gammaj.elem( valid_elems ) ), gamma );

                if( gammatilde < gamma )
                {
                    gamma = gammatilde;
                    drop = find( gammaj == gammatilde );
                }
            }
        }

        beta.elem( active.get_active( ) ) = beta.elem( active.get_active( ) ) + w * gamma;
        mu = mu + gamma * u;
        double model_var = sum( pow( phenotype - mu, 2 ) ) / (n - 1 - active.size( ) );

        if( lasso && drop.n_elem > 0 )
        {
            beta.elem( active.get_active( ).elem( drop ) ).fill( 0.0 );
            
            active.drop( drop[ 0 ] );

            /* Don't report drops */
            continue;
        }
        
        /**
         * Compute value of the knot.
         */
        arma::vec r = phenotype - mu;
        double lambda = arma::max( s % ( X_active.t( ) * r ) );

        /* Compute p-values and add path of previous coefficients */
        arma::uvec cur_active = active.get_active( );
        float tsum = sum( abs( beta ) );
        for(int j = 0; j < cur_active.n_elem - 1; j++)
        {
            if( cur_active[ j ] != max_index )
            {
                path.push_back( lars_path( i + 1,
                                           variable_set.get_name( cur_active[ j ] ),
                                           beta[ cur_active[ j ] ],
                                           0.0,
                                           1.0,
                                           lambda,
                                           tsum,
                                           1.0 - model_var / pheno_var ) );
            }
        }
        
        /* Calculate p for new beta and add last component of the path */
        double cur_cor = dot( mu, phenotype );

        /* Compute h0 */
        double prev_cor = base_cor;
        if( prev_active.n_elem > 0 )
        {
            arma::mat X_h0 = variable_set.get_active( prev_active );
            arma::vec beta_h0 = optimize_lars_gd( X_h0, phenotype, lambda, beta.elem( prev_active ) );
            prev_cor = dot( X_h0 * beta_h0, phenotype );
        }

        double T = (cur_cor - prev_cor ) / model_var;
        double p = 0.9999;
        if( T > 0.0 )
        {
            p = 1 - exp_cdf( T, 1.0 );
        }

        path.push_back( lars_path( i + 1,
                        variable_set.get_name( max_index ),
                        beta[ max_index ],
                        T,
                        p,
                        lambda,
                        tsum,
                        1.0 - model_var / pheno_var ) );

        last_p = p;
        i++;
    }

    return path;
}

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG ); 
 
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "-m", "--max-variables" ).help( "Maximum number of variables in the model" ).set_default( 10 );
    parser.add_option( "-t", "--threshold" ).help( "Stop after a variable has a p-value less than this threshold." ).set_default( 1.0 );
    parser.add_option( "--only-pvalues" ).help( "Only output the beta that enters in each step along with its p-value." ).action( "store_true" );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 1 )
    {
        parser.print_help( );
        exit( 1 );
    }

    std::ios_base::sync_with_stdio( false );
    
    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( parser.args( )[ 0 ] );
    genotype_matrix_ptr genotypes = create_genotype_matrix( genotype_file );
    std::vector<std::string> order = genotype_file->get_sample_iids( );

    /* Make error streams separate from stdout */
    arma::set_stream_err1( std::cerr );
    arma::set_stream_err2( std::cerr );

    /* Parse phenotypes */
    arma::uvec cov_missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::uvec pheno_missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::vec phenotype;
    arma::mat cov;
    std::vector<std::string> cov_names;
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        phenotype = parse_phenotypes( phenotype_file, pheno_missing, order, options[ "mpheno" ] );
    }
    else
    {
        phenotype = create_phenotype_vector( genotype_file->get_samples( ), pheno_missing );
    }
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        cov = parse_covariate_matrix( covariate_file, cov_missing, order, &cov_names );
    }

    /* Open output stream */
    std::ofstream output_file;
    if( options.is_set( "out" ) )
    {
        output_file.open( options[ "out" ].c_str( ) );
    }
    std::ostream &out = options.is_set( "out" ) ? output_file : std::cout;

    bool only_pvalues = options.is_set( "only_pvalues" );

    gene_environment variable_set( genotypes, cov, phenotype, cov_names );
    std::vector<lars_path> path = lars( variable_set, (int) options.get( "max_variables" ), (double) options.get( "threshold" ) );
    out << "step\tvariable\tbeta\tT\tp\tlambda\tbeta_sum\texplained_var\n";
    for(int i = 0; i < path.size( ); i++)
    {
        lars_path step = path[ i ];
        if( !only_pvalues || (only_pvalues && step.m_pvalue < 1.0) )
        {
            out << step.m_pass << "\t" <<
                step.m_variable << "\t" <<
                step.m_beta << "\t" <<
                step.m_T << "\t" <<
                step.m_pvalue << "\t" <<
                step.m_lambda <<  "\t" <<
                step.m_beta_sum <<  "\t" <<
                step.m_explained_var << "\n";
        }
    }

    return 0;
}
