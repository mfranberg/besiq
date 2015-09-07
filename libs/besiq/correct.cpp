#include <algorithm>
#include <iostream>

#include <glm/models/binomial.hpp>
#include <besiq/io/resultfile.hpp>
#include <besiq/io/metaresult.hpp>
#include <besiq/method/scaleinv_method.hpp>

#include <besiq/correct.hpp>

struct heap_result
{
    std::pair<std::string, std::string> variant_pair;
    float pvalue;

    bool operator<(const heap_result &other)
    {
        return other.pvalue > pvalue;
    }
};

void
run_top(metaresultfile *result, float alpha, uint64_t num_top, size_t column, const std::string &output_path)
{
    std::ostream &output = std::cout;
    std::vector<std::string> header = result->get_header( );
    output << "snp1 snp2\tP\n";

    std::vector<heap_result> heap;
    std::make_heap( heap.begin( ), heap.end( ) );

    std::pair<std::string, std::string> pair;
    float *values = new float[ header.size( ) ];
    float cur_min = FLT_MAX;
    while( result->read( &pair, values ) )
    {
        heap_result res;
        res.pvalue = values[ column ];
        if( res.pvalue == result_get_missing( ) || res.pvalue > cur_min )
        {
            continue;
        }

        res.variant_pair = pair;
        heap.push_back( res );
        std::push_heap( heap.begin( ), heap.end( ) );
        if( heap.size( ) > num_top )
        {
            std::pop_heap( heap.begin( ), heap.end( ) );
            float new_min = heap.back( ).pvalue;
            heap.pop_back( );

            cur_min = std::min( cur_min, new_min );
        }
    }
    std::sort_heap( heap.begin( ), heap.end( ) );
    for(int i = 0; i < heap.size( ); i++)
    {
        output << heap[ i ].variant_pair.first << " " << heap[ i ].variant_pair.second;
        output << "\t" << heap[ i ].pvalue << "\n";
    }
}

void
run_bonferroni(metaresultfile *result, float alpha, uint64_t num_tests, size_t column, const std::string &output_path)
{
    if( num_tests == 0 )
    {
        num_tests = result->num_pairs( );
    }

    std::ostream &output = std::cout;
    std::vector<std::string> header = result->get_header( );
    output << "snp1 snp2";
    for(int i = 0; i < header.size( ); i++)
    {
        output << "\t" << header[ i ];
    }
    output << "\tP_adjusted\n";

    std::pair<std::string, std::string> pair;
    float *values = new float[ header.size( ) ];
    while( result->read( &pair, values ) )
    {
        float p = values[ column ];
        if( p == result_get_missing( ) )
        {
            continue;
        }

        float adjusted_p = std::min( p * num_tests, 1.0f );

        if( adjusted_p <= alpha )
        {
            output << pair.first << " " << pair.second;
            for(int i = 0; i < header.size( ); i++)
            {
                output << "\t" << values[ i ];
            }
            output << "\t" << adjusted_p << "\n";
        }
    }
}

resultfile *
do_common_stages(metaresultfile *result, const std::vector<std::string> &snp_names, const correction_options &options, const std::string &output_path)
{
    float *values = new float[ result->get_header( ).size( ) ];
    char const *levels[] = { "1", "2", "3", "4" };
    std::string filename = output_path + std::string( ".level" ) + levels[ 0 ];
    bresultfile *stage_file = new bresultfile( filename, snp_names );
    if( stage_file == NULL || !stage_file->open( ) )
    {
        return NULL;
    }

    stage_file->set_header( result->get_header( ) );

    std::vector<uint64_t> num_tests( options.num_tests );
    if( num_tests[ 0 ] == 0 )
    {
        num_tests[ 0 ] = result->num_pairs( );
    }
    std::pair<std::string, std::string> pair;
    while( result->read( &pair, values ) )
    {
        if( values[ 0 ] == result_get_missing( ) )
        {
            continue;
        }
        
        float adjusted = std::min( 1.0f, values[ 0 ] * num_tests[ 0 ] / options.weight[ 0 ] );
        if( adjusted <= options.alpha )
        {
            values[ 0 ] = adjusted;
            stage_file->write( pair, values );
        }
    }
    delete stage_file;

    for(int i = 1; i < 3; i++)
    {
        bresultfile *prev_file = new bresultfile( filename );
        if( prev_file == NULL || !prev_file->open( ) )
        {
            return NULL;
        }

        filename = output_path + std::string( ".level" ) + levels[ i ];
        stage_file = new bresultfile( filename, snp_names );
        if( stage_file == NULL || !stage_file->open( ) )
        {
            return NULL;
        }

        stage_file->set_header( result->get_header( ) );
        
        if( num_tests[ i ] == 0 )
        {
            num_tests[ i ] = prev_file->num_pairs( );
        }

        while( prev_file->read( &pair, values ) )
        {
            if( values[ i ] == result_get_missing( ) )
            {
                continue;
            }
            
            float adjusted = std::min( values[ i ] * num_tests[ i ] / options.weight[ i ], 1.0f );
            if( adjusted <= options.alpha )
            {
                values[ i ] = adjusted;
                stage_file->write( pair, values );
            }
        }

        delete stage_file;
        delete prev_file;
    }

    delete values;

    bresultfile *last_stage = new bresultfile( filename );
    if( last_stage != NULL && last_stage->open( ) )
    {
        return last_stage;
    }
    else
    {
        return NULL;
    }
}

void
do_last_stage(resultfile *last_stage, const correction_options &options, genotype_matrix_ptr genotypes, method_data_ptr data, const std::string &output_path)
{
    std::ostream &output = std::cout;
    std::vector<std::string> header = last_stage->get_header( );

    model_matrix *model_matrix = new factor_matrix( data->covariate_matrix, data->phenotype.n_elem );
    scaleinv_method method( data, *model_matrix, options.model == "normal" );
    std::vector<std::string> method_header = method.init( ); 

    output << "snp1 snp2";
    for(int i = 0; i < method_header.size( ); i++)
    {
        output << "\tP_" << method_header[ i ];
    }
    output << "\tP_combined\n";

    uint64_t num_tests = options.num_tests[ 3 ];
    if( num_tests == 0 )
    {
        num_tests = last_stage->num_pairs( );
    }

    std::pair<std::string, std::string> pair;
    float *values = new float[ header.size( ) ];
    float *method_values = new float[ method_header.size( ) ];
    while( last_stage->read( &pair, values ) )
    {
        std::fill( method_values, method_values + method_header.size( ), result_get_missing( ) );
        /* Skip N and last p-value */
        float pre_p = *std::max_element( values, values + header.size( ) - 2 );
        float min_p = 1.0;
        float max_p = 0.0;
        
        snp_row const *snp1 = genotypes->get_row( pair.first );
        snp_row const *snp2 = genotypes->get_row( pair.second );
        method.run( *snp1, *snp2, method_values );

        std::vector<float> p_values;
        for(int i = 0; i < method_header.size( ); i++)
        {
            float adjusted_p = result_get_missing( );    
            if( method_values[ i ] != result_get_missing( ) )
            {
                adjusted_p = std::min( std::max( method_values[ i ] * num_tests / options.weight[ header.size( ) - 2 ], pre_p ), 1.0f );
                min_p = std::min( min_p, adjusted_p );
                max_p = std::max( max_p, adjusted_p );
            }

            p_values.push_back( adjusted_p );
        }

        if( min_p > options.alpha )
        {
            continue;
        }

        output << pair.first << " " << pair.second;
        bool any_missing = false;
        for(int i = 0; i < p_values.size( ); i++)
        {
            if( p_values[ i ] != result_get_missing( ) )
            {
                output << "\t" << p_values[ i ];
            }
            else
            {
                any_missing = true;
                output << "\tNA";
            }
        }

        if( !any_missing )
        {
            output << "\t" << max_p << "\n";
        }
        else
        {
            output << "\tNA\n";
        }
    }
    
    delete model_matrix;
    delete values;
    delete method_values;
}

void
run_static(metaresultfile *result, genotype_matrix_ptr genotypes, method_data_ptr data, const correction_options &options, const std::string &output_path)
{
    resultfile *last_stage = do_common_stages( result, genotypes->get_snp_names( ), options, output_path );
    if( last_stage == NULL )
    {
        return;
    }

    do_last_stage( last_stage, options, genotypes, data, output_path );

    delete last_stage;
}

void
run_adaptive(metaresultfile *result, genotype_matrix_ptr genotypes, method_data_ptr data, const correction_options &options, const std::string &output_path)
{
    correction_options adaptive_options( options );
    for(int i = 0; i < adaptive_options.num_tests.size( ); i++)
    {
        adaptive_options.num_tests[ i ] = 0;
    }

    run_static( result, genotypes, data, adaptive_options, output_path );
}

