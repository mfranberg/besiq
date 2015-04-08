#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <algorithm>
#include <string>

extern char **environ;

struct command
{
    const char *name;
    const char *description;
};

struct command g_commands[] =
{
    { "pairs", "Generate pairs for a gene-gene analysis." },
    { "view", "Display a binary result file" },
    { "glm", "Run a GLM model possibly with covariates." },
    { "scaleinv", "Scale-invariant tests of interaction." },
    { "stagewise", "Run a stage-wise method." },
    { "wald", "Perform a 'fast' wald test with main effects." },
    { "bayes", "Run a stage-wise method." },
    { "loglinear", "Run the log-linear method without main effects." },
    { "caseonly", "Run the case only test." },
    { "correct", "Multiple testing correction." },
    { "env", "Run a gene-environment test." },
    { NULL, NULL }
};

void print_help()
{
    printf( "Usage: bayesic COMMAND [ARGS]...\n\n" );
    printf( "  A collection of tools for inferring genetic interactions\n" );
    printf( "  that include generation of pairs to test, different\n" );
    printf( "  tests for GxG and GxE interaction, multiple testing\n" );
    printf( "  correction, and viewing binary result files.\n\n" );
    printf( "Commands:\n" );
    
    size_t max_len = 0;
    for(int i = 0; g_commands[ i ].name != NULL; i++)
    {
        if( strlen( g_commands[ i ].name ) > max_len ) 
        {
            max_len = strlen( g_commands[ i ].name );
        }
    }
    
    for(int i = 0; g_commands[ i ].name != NULL; i++)
    {
        printf( "  %s", g_commands[ i ].name );
        for(int j = 0; j < max_len - strlen( g_commands[ i ].name ) + 2; j++ )
        {
            putchar( ' ' );
        };
        printf( "%s\n", g_commands[ i ].description );
    }
}

bool
find_command(const char *cmd)
{
    for(int i = 0; g_commands[ i ].name != NULL; i++)
    {
        if( strcmp( g_commands[ i ].name, cmd ) == 0 )
        {
            return true;
        }
    }

    return false;
}

std::string dirname(std::string source)
{
    if( source.size( ) <= 1 )
    {
        return source;
    }

    source.erase( std::find( source.rbegin( ), source.rend( ), '/' ).base( ), source.end( ) );

    return source;
}

int main(int argc, char *argv[])
{
    if( argc < 2 || !find_command( argv[ 1 ] ) )
    {
        print_help( );
        exit( 1 );
    }

    if( !find_command( argv[ 1 ] ) )
    {
        fprintf( stderr, "bayesic: error: No such command %s.\n", argv[ 1 ] );
        exit( 1 );
    }

    std::string cmd_path( argv[ 0 ] );
    std::string cmd_dir = dirname( cmd_path );
    std::string new_cmd_path = cmd_dir + "bayesic-" + std::string( argv[ 1 ] );

    argv[ 1 ] = strdup( new_cmd_path.c_str( ) );
    if( execve( new_cmd_path.c_str( ), &argv[ 1 ], environ ) == -1 )
    {
        fprintf( stderr, "bayesic: error: Could not find executable %s.\n", new_cmd_path.c_str( ) );
        exit( 1 );
    }
}
