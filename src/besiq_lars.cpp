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

struct lars_path
{
    lars_path(int p, std::string &v, double b, double pval, double t, double e)
        : pass( p ),
          variable( v ),
          beta( b ),
          pvalue( pval ),
          tsum( t ),
          exp_var( e )
    {
    }

    int pass;
    std::string variable;
    double beta;
    double pvalue;
    double tsum;
    double exp_var;
};

void
fill_missing_genotypes(genotype_matrix_ptr genotypes)
{
    #pragma omp parallel for
    for(int i = 0; i < genotypes->size( ); i++)
    {
        snp_row &row = genotypes->get_row( i );
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
fill_missing_phenotypes(arma::vec &phenotypes)
{
    double mu = 0.0;
    int total = 0;
    for(int i = 0; i < phenotypes.n_elem; i++)
    {
        if( phenotypes[ i ] == phenotypes[ i ] )
        {
            mu += phenotypes[ i ];
            total++;
        }
    }
    mu = mu / total;

    for(int i = 0; i < phenotypes.n_elem; i++)
    {
        if( phenotypes[ i ] != phenotypes[ i ] )
        {
            phenotypes[ i ] = mu;
        }
    }
}

void
fill_missing_cov(arma::mat &cov)
{
    for(int i = 0; i < cov.n_cols; i++)
    {
        double mu = 0.0;
        int total = 0;
        for(int j = 0; j < cov.n_rows; j++)
        {
            // Check NaN, may fail in some compilers
            if( cov( j, i ) == cov( j, i ) )
            {
                mu += cov( j, i );
                total++;
            }
        }
        mu = mu / total;

        for(int j = 0; j < cov.n_rows; j++)
        {
            // Check NaN, may fail in some compilers
            if( cov( j, i ) != cov( j, i ) )
            {
                cov( j, i ) = mu;
            }
        }
    }
}

arma::vec
nice_division(arma::vec &a, arma::vec &b)
{
    b.elem( find( b < DBL_MIN ) ).fill( DBL_MIN );
    arma::vec x = a / b;
    x.elem( find( x <= 0 ) ).fill( DBL_MAX );
    
    return x;
}

void
calculate_cor(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &phenotype, const arma::vec &mu, const arma::vec &mean, const arma::vec &sd, arma::vec &c)
{
    int var = 0;
    arma::vec adjusted_pheno = phenotype - mu;
    #pragma omp parallel for
    for(int i = 0; i < genotypes->size( ); i++)
    {
        const snp_row &row = genotypes->get_row( i );
        double cor = 0.0;
        double m = mean[ i ];
        double s = sd[ i ];
        for(int j = 0; j < row.size( ); j++)
        {
            cor += ( ( row[ j ] - m ) / s ) * adjusted_pheno[ j ];
        }
        c[ i ] = cor;
    }

    var = genotypes->size( );
    for(int j = 0; j < cov.n_cols; j++)
    {
        double cor = dot( ( cov.col( j ) - mean[ var + j ] ) / sd[ var + j ],  adjusted_pheno );
        c[ var + j ] = cor;
    }
}

std::vector<std::string>
create_names(const std::vector<std::string> &locus_names, const std::vector<std::string> &cov_names)
{
    std::vector<std::string> names( locus_names );
    if( cov_names.size( ) > 0 )
    {
        names.insert( names.end( ), cov_names.begin( ) + 2, cov_names.end( ) );
    }

    return names;
}
arma::vec eig_prod(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &u, const arma::vec &mean, const arma::vec &sd)
{
    arma::vec a = arma::zeros<arma::vec>( genotypes->size( ) + cov.n_cols );

    size_t var = 0;
    #pragma omp parallel for
    for(int i = 0; i < genotypes->size( ); i++)
    {
        const snp_row &row = genotypes->get_row( i );
        double m = mean[ i ];
        double s = sd[ i ];
        for(int j = 0; j < row.size( ); j++)
        {
            a[ i ] += ( ( row[ j ] - m ) / s ) * u[ j ];
        }
    }

    var = genotypes->size( );
    for(int i = 0; i < cov.n_cols; i++)
    {
        a[ var + i ] = dot( ( cov.col( i ) - mean[ var + i ] ) / sd[ var + i ], u );
    }

    return a;
}

arma::mat get_active(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::uvec &active, size_t n, const arma::vec &mean, const arma::vec &sd)
{
    arma::mat X( n, active.n_elem );
    size_t var = 0;
    for(int i = 0; i < active.n_elem; i++)
    {
        double m = mean[ active[ i ] ];
        double s = sd[ active[ i ] ];
        if( active[ i ] < genotypes->size( ) )
        {
            const snp_row &row = genotypes->get_row( active[ i ] );
            for(int j = 0; j < row.size( ); j++)
            {
                X( j, var ) = ( row[ j ] - m ) / s;
            }
        }
        else
        {
            X.col( var ) = ( cov.col( active[ i ] - genotypes->size( ) ) - m ) / s;
        }
        var++;
    }

    return X;
}

void compute_mean_sd(genotype_matrix_ptr genotypes, const arma::mat &cov, arma::vec &mean, arma::vec &sd, size_t n)
{
    size_t var = 0;
    for(int i = 0; i < genotypes->size( ); i++)
    {
        snp_row &row = genotypes->get_row( i );
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

        mean[ var ] = m;
        sd[ var ] = sqrt( s );

        var++;
    }

    for(int i = 0; i < cov.n_cols; i++)
    {
        double m = arma::mean( cov.col( i ) );
        double s = arma::accu( pow( cov.col( i ) - m, 2 ) );

        mean[ var ] = m;
        sd[ var ] = sqrt( s );

        var++;
    }
}

arma::uvec get_complement(const arma::uvec &active, size_t m)
{
    arma::uvec all = arma::zeros<arma::uvec>( m );
    all.elem( active ).fill( 1 );
    return find( all == 0 );
}

std::vector<lars_path>
lars(genotype_matrix_ptr genotypes, arma::mat &cov, arma::vec &phenotype, const std::vector<std::string> &cov_names, size_t max_vars = 15, bool lasso = true, bool only_p = false)
{
    fill_missing_genotypes( genotypes );
    fill_missing_cov( cov );
    fill_missing_phenotypes( phenotype );

    phenotype = phenotype - mean( phenotype );
    size_t n = phenotype.n_elem;
    size_t m = genotypes->get_snp_names( ).size( ) + cov.n_cols;
    arma::vec mu = arma::zeros<arma::vec>( n );
    arma::vec beta = arma::zeros<arma::vec>( m );
    std::vector<lars_path> path;
    arma::vec c = arma::zeros<arma::vec>( m );
    std::vector<std::string> names = create_names( genotypes->get_snp_names( ), cov_names );

    arma::vec mean = arma::zeros<arma::vec>( m );
    arma::vec sd = arma::zeros<arma::vec>( m );
    compute_mean_sd( genotypes, cov, mean, sd, n );

    double pheno_mean = sum( phenotype ) / n;
    double pheno_var = sum( pow( phenotype - pheno_mean, 2 ) ) / (n - 1);
    double prev_cor = sum( phenotype * pheno_mean );

    arma::uvec drop;
    arma::uvec dropinv;
    arma::uvec active;
    arma::uvec inactive = arma::linspace<arma::uvec>( 0, m - 1, m );
    double eps = 1e-10;


    int i = 0;
    while(active.n_elem < std::min( max_vars, m ) )
    {
        calculate_cor( genotypes, cov, phenotype, mu, mean, sd, c );
        arma::vec cabs = abs( c );
        double C = max( cabs.elem( inactive ) );
        
        if( C < eps )
        {
            break;
        }
        
        if( drop.n_elem == 0 )
        {
            arma::uvec new_vars = find( cabs.elem( inactive ) >= C - eps );
            active.resize( active.n_elem + new_vars.n_elem );
            active.tail( new_vars.n_elem ) = inactive.elem( new_vars );
            inactive = get_complement( active, m );
        }

        /* Create the active matrix */
        arma::mat X_active = get_active( genotypes, cov, active, n, mean, sd );

        /* Equation 2.4 */
        arma::vec s = sign( c.elem( active ) );

        /* Equation 2.5 */
        arma::mat G = X_active.t( ) * X_active;
        arma::mat Ginv = inv( G );
        double A = 1.0 / sqrt( arma::accu( s.t( ) * Ginv * s ) );

        /* Equation 2.6 */
        arma::vec w = A * Ginv * s;
        arma::vec u = X_active * w;
        
        double gamma = C / A;
        if( active.n_elem < m )
        { 
            /* Equation 2.11 */
            arma::vec a = eig_prod( genotypes, cov, u, mean, sd );

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
            arma::vec gammaj = -beta.elem( active ) / w;

            /* Equation 3.5 */
            arma::uvec valid_elems = arma::find( gammaj > eps );
            if( valid_elems.n_elem > 0 )
            {
                double gammatilde = std::min( arma::min( gammaj.elem( valid_elems ) ), gamma );

                if( gammatilde < gamma )
                {
                    gamma = gammatilde;
                    drop = find( gammaj == gammatilde );
                    dropinv = find( gammaj != gammatilde );
                }
            }
        }

        beta.elem( active ) = beta.elem( active ) + w * gamma;
        mu = mu + gamma * u;
        double model_var = sum( pow( phenotype - mu, 2 ) ) / (n - 1 - active.n_elem);

        if( lasso && drop.n_elem > 0 )
        {
            beta.elem( active.elem( drop ) ).fill( 0.0 );
            active = active.elem( dropinv );
            inactive = get_complement( active, m );
        }

        float tsum = sum( abs( beta ) );
        for(int j = 0; j < active.n_elem - 1; j++)
        {
            path.push_back( lars_path( i + 1, names[ active[ j ] ], beta[ active[ j ] ], 1.0, tsum, 1.0 - model_var / pheno_var ) );
        }
        
        /* Calculate p for new beta and add last component of the path */
        double cur_cor = dot( mu, phenotype );
        double T = (cur_cor - prev_cor ) / pheno_var;
        double p = 1 - exp_cdf( T, 1 );
        prev_cor = cur_cor;
        path.push_back( lars_path( i + 1, names[ active[ active.n_elem - 1 ] ], beta[ active[ active.n_elem - 1 ] ], p, tsum, 1.0 - model_var / pheno_var ) );

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

    bool only_pvalues = options.is_set( "only_pvalues" );
    std::vector<lars_path> path = lars( genotypes, cov, phenotype, cov_names, (int) options.get( "max_variables" ) );
    std::cout << "step\tvariable\tbeta\tp\tsum\texplained_var\n";
    for(int i = 0; i < path.size( ); i++)
    {
        lars_path step = path[ i ];
        if( !only_pvalues || (only_pvalues && step.pvalue < 1.0) )
        {
            std::cout << step.pass << "\t" << step.variable << "\t" << step.beta << "\t" << step.pvalue << "\t" << step.tsum <<  "\t" << step.exp_var << "\n";
        }
    }

    return 0;
}
