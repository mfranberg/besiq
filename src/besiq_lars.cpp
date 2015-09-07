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
    lars_path(int p, std::string &v, float b, float t)
        : pass( p ),
          variable( v ),
          beta( b ),
          tsum( t )
    {
    }

    int pass;
    std::string variable;
    float beta;
    float tsum;
};

void
fill_missing_genotypes(genotype_matrix_ptr genotypes)
{
    for(int i = 0; i < genotypes->size( ); i++)
    {
        snp_row &row = genotypes->get_row( i );
        double maf = compute_real_maf( row );
        int mu = round( 2*maf*(1-maf) + 2*maf*maf );
        for(int j = 0; j < row.size( ); j++)
        {
            if( row[ j ] == 3 )
            {
                row.assign( j, mu );
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
        if( phenotypes[ i ] != arma::datum::nan )
        {
            mu += phenotypes[ i ];
            total++;
        }
    }
    mu = mu / total;

    for(int i = 0; i < phenotypes.n_elem; i++)
    {
        if( phenotypes[ i ] == arma::datum::nan )
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
        for(int j = 0; j < cov.n_rows; i++)
        {
            if( cov( i, j ) != arma::datum::nan )
            {
                mu += cov( i, j );
                total++;
            }
        }
        mu = mu / total;

        for(int j = 0; j < cov.n_rows; j++)
        {
            if( cov( i, j ) == arma::datum::nan )
            {
                cov( i, j ) = mu;
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
    for(int i = 0; i < genotypes->size( ); i++)
    {
        const snp_row &row = genotypes->get_row( i );
        double cor = 0.0;
        double m = mean[ var ];
        double s = sd[ var ];
        for(int j = 0; j < row.size( ); j++)
        {
            cor += ( ( row[ j ] - m ) / s ) * adjusted_pheno[ j ];
        }
        c[ var ] = cor;
        var++;
    }

    for(int j = 0; j < cov.n_cols; j++)
    {
        double cor = dot( ( cov.col( j ) - mean[ var ] ) / sd[ var ],  adjusted_pheno );
        c[ var ] = cor;
        var++;
    }
}

std::vector<std::string>
create_names(const std::vector<std::string> &locus_names, const std::vector<std::string> &cov_names)
{
    std::vector<std::string> names( locus_names );
    names.insert( names.end( ), cov_names.begin( ), cov_names.end( )) ;

    return names;
}
arma::vec eig_prod(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &u, const arma::vec &mean, const arma::vec &sd)
{
    arma::vec a = arma::zeros<arma::vec>( genotypes->size( ) );

    size_t var = 0;
    for(int i = 0; i < genotypes->size( ) + cov.n_cols; i++)
    {
        const snp_row &row = genotypes->get_row( i );
        double m = mean[ var ];
        double s = sd[ var ];
        for(int j = 0; j < row.size( ); j++)
        {
            a[ i ] += ( ( row[ j ] - m ) / s ) * u[ j ];
        }
        var++;
    }
    for(int i = 0; i < cov.n_cols; i++)
    {
        a[ i ] = dot( ( cov.col( i ) - mean[ var ] ) / sd[ var ], u );
        var++;
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
        sd[ var ] = sqrt( s / n );

        var++;
    }

    for(int i = 0; i < cov.n_cols; i++)
    {
        double m = arma::mean( cov.col( i ) );
        double s = arma::accu( pow( cov.col( i ) - m, 2 ) );

        mean[ var ] = m;
        sd[ var ] = sqrt( s / n );

        var++;
    }
}

std::vector<lars_path>
lars(genotype_matrix_ptr genotypes, arma::mat &cov, arma::vec &phenotype, const std::vector<std::string> &cov_names, int max_vars = 15)
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
    std::set<int> active;
    arma::vec c = arma::zeros<arma::vec>( m );
    std::vector<std::string> names = create_names( genotypes->get_snp_names( ), cov_names );

    arma::vec mean = arma::zeros<arma::vec>( m );
    arma::vec sd = arma::zeros<arma::vec>( m );
    compute_mean_sd( genotypes, cov, mean, sd, n );

    arma::uvec drop;

    for(int i = 0; i < max_vars; i++)
    {
        calculate_cor( genotypes, cov, phenotype, mu, mean, sd, c );
        arma::vec cabs = abs( c );
        double C = max( cabs );
        
        if( C == 0 || C < 1e-10 )
        {
            break;
        }
        
        uvec active = find( abs( cabs - C ) < 1e-10 + 1e-10 * C );
        uvec inactive = find( abs( cabs - C ) > 1e-10 + 1e-10 * C );
        if( active.n_elem == m )
        { 
            break;
        }

        /* Create the active matrix */
        arma::mat X = get_active( genotypes, cov, active, n, mean, sd );

        arma::vec s = sign( c.elem( active ) );
        arma::mat X_active = X * diagmat( s );

        arma::mat G = X_active.t( ) * X_active;
        arma::mat Ginv = inv( G );

        double A = sqrt( arma::accu( Ginv ) );

        arma::vec w = A * Ginv * arma::ones<arma::vec>( active.n_elem );
        arma::vec u = X_active * w;
        
        arma::vec a = eig_prod( genotypes, cov, u, mean, sd );

        arma::vec cc = c.elem( inactive );
        arma::vec ac = a.elem( inactive );

        arma::vec a1 = C - cc;
        arma::vec b1 = A - ac;
        arma::vec gn = nice_division( a1, b1 );

        arma::vec a2 = C + cc;
        arma::vec b2 = A + ac;
        arma::vec gp = nice_division( a2, b2 );

        double gamma = std::min( gn.min( ), gp.min( ) );
        
        /*if( lasso )
        {
            double gammaj = -beta.elem( active ) / w;
            double gammatilde = std::min( min( gammaj.elem( find( gammaj > 1e-11 ) ) ), gamma );

            if( gammatilde < gamma )
            {
                gamma = gammatilde;
                drop = find( gammaj = gammatilde );
            }
        }*/

        beta.elem( active ) = beta.elem( active ) + s % w * gamma;
        mu = mu + gamma * u;
        float tsum = sum( abs( beta ) );
        for(int j = 0; j < active.n_elem; j++)
        {
            path.push_back( lars_path( i, names[ active[ j ] ], beta[ active[ j ] ], tsum ) );
        }
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

    std::vector<lars_path> path = lars( genotypes, cov, phenotype, cov_names, (int) options.get( "max_variables" ) );
    std::cout << "step\tvariable\tbeta\ttsum\n";
    for(int i = 0; i < path.size( ); i++)
    {
        lars_path step = path[ i ];
        std::cout << step.pass << "\t" << step.variable << "\t" << step.beta << "\t" << step.tsum << "\n";
    }

    return 0;
}
