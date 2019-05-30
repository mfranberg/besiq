#include <iostream>

#include <armadillo>

#include <glm/models/binomial.hpp>
#include <glm/models/normal.hpp>
#include <cpp-argparse/OptionParser.h>

#include <besiq/method/scaleinv_method.hpp>
#include <besiq/io/covariates.hpp>
#include <besiq/io/resultfile.hpp>
#include <besiq/stats/snp_count.hpp>
#include <plink/imputed.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-imputed [OPTIONS] gen1 gen2 pheno";
const std::string DESCRIPTION = "Computes p-values for imputed data";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

ssize_t find_variant_index(const std::vector<imputed_info> &info, const std::string &name)
{
    for(int i = 0; i < info.size( ); i++)
    {
        if( name == info[ i ].name )
        {
            return i;
        }
    }
    return -1;
}

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG ); 
 
    char const* const model_choices[] = { "binomial", "normal" };

    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );
    parser.add_option( "-w", "--window" ).help( "Number of variants around the best pair to look at." ).set_default( 10 );
    parser.add_option( "-i", "--info" ).help( "Info threshold" ).set_default( 0.8 );
    parser.add_option( "-m", "--maf" ).help( "MAF threshold" ).set_default( 0.05 );
    parser.add_option( "--snp1" ).help( "The first variant." );
    parser.add_option( "--snp2" ).help( "The second variant." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );
    parser.add_option( "--call" ).help( "Call rate for hard calls (default = 0.8)." ).set_default( 0.8 );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 3 )
    {
        parser.print_help( );
        exit( 1 );
    }

    std::ios_base::sync_with_stdio( false );
    float call_rate = (float) options.get( "call" );
    int window_size = (int) options.get( "window" );
    float info_threshold = (float) options.get( "info" );
    float maf_threshold = (float) options.get( "maf" );

    method_data_ptr data( new method_data( ) );
    imputed_data imputed1 = parse_imputed_data( parser.args( )[ 0 ], call_rate );
    imputed_data imputed2 = parse_imputed_data( parser.args( )[ 1 ], call_rate );
    data->missing = arma::zeros<arma::uvec>( imputed1.samples.size( ) );
    
    int variant1 = find_variant_index( imputed1.info, options[ "snp1" ] );
    int variant2 = find_variant_index( imputed2.info, options[ "snp2" ] );
    if( variant1 == -1 || variant2 == -1 )
    {
        std::cerr << "besiq-imputed: error: Could not find the specified variants: '" << options[ "snp1" ] << "' or '" << options[ "snp2" ] << "'\n";
        exit( 1 );
    }

    std::ifstream pheno_file( parser.args( )[ 2 ] );
    data->phenotype = parse_phenotypes( pheno_file, data->missing, imputed1.samples, options[ "mpheno" ] );
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        data->covariate_matrix = parse_covariate_matrix( covariate_file, data->missing, imputed1.samples );
    }
    
    /* Open results. */
    resultfile *result_file = NULL;
    if( options.is_set( "out" ) )
    {
        result_file = new tresultfile( options[ "out" ], "w" );
    }
    else
    {
        result_file = new tresultfile( "-", "w" );
    }
    if( result_file == NULL || !result_file->open( ) )
    {
        std::cerr << "besiq: error: Can not open result file." << std::endl;
        exit( 1 );
    }
    
    model_matrix *model_matrix = make_model_matrix( "factor", data->covariate_matrix, data->phenotype.n_elem );
    method_type *m = NULL;
    if( options[ "model" ] == "binomial" )
    {
        m = new scaleinv_method( data, *model_matrix, false );
    }
    else if( options[ "model" ] == "normal" )
    {
        m = new scaleinv_method( data, *model_matrix, true );
    }

    std::vector<std::string> method_header = m->init( );
    size_t method_size = method_header.size( );
    method_header.push_back( "maf1" );
    method_header.push_back( "maf2" );
    method_header.push_back( "info1" );
    method_header.push_back( "info2" );
    method_header.push_back( "N" );
    result_file->set_header( method_header );

    float *output = new float[ method_header.size( ) ];

    int v1_start = std::max( variant1 - window_size, 0 );
    int v1_end = std::min( variant1 + window_size, (int) imputed1.genotypes->size( ) );
    int v2_start = std::max( variant2 - window_size, 0 );
    int v2_end = std::min( variant2 + window_size, (int) imputed2.genotypes->size( ) );

    for(int i = v1_start; i < v1_end; i++)
    {
        for(int j = v2_start; j < v2_end; j++)
        {
            if( imputed1.info[ i ].name == imputed2.info[ j ].name )
            {
                continue;
            }

            std::fill( output, output + method_header.size( ), result_get_missing( ) );
        

            float info1 = imputed1.info[ i ].info;
            float info2 = imputed2.info[ j ].info;

            const snp_row &row1 = (*imputed1.genotypes)[ i ];
            const snp_row &row2 = (*imputed2.genotypes)[ j ];
            
            float maf1 = compute_real_maf( row1 );
            float maf2 = compute_real_maf( row2 );

            if( maf1 > maf_threshold && maf2 > maf_threshold && info1 > info_threshold && info2 > info_threshold )
            {
                m->run( row1, row2, output );
            
                output[ method_size ] = maf1;
                output[ method_size + 1 ] = maf2;
                output[ method_size + 2 ] = info1;
                output[ method_size + 3 ] = info2;
                output[ method_size + 4 ] = m->num_ok_samples( row1, row2 );

                std::pair<std::string, std::string> pair( imputed1.info[ i ].name, imputed2.info[ j ].name );
                result_file->write( pair, output );
            }
        }
    }

    return 0;
}
