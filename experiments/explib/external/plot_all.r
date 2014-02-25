library( ggplot2 )
library( gtable )

argv = commandArgs( trailingOnly = TRUE )
if ( length( argv ) != 5 )
{
    message( "error: Wrong number of arguments." )
    message( "Usage: plot_model_power xlabel ylabel title model_input output_file" )
    quit( )
}

xlabel = argv[ 1 ]
ylabel = argv[ 2 ]
title = argv[ 3 ]

model_power = read.table( argv[ 4 ], header = FALSE )
names( model_power ) = c( "method", "model", "heritability", "power", "lower", "upper" )

output_file = argv[ 5 ]

pdf( output_file, width = 2 * 6.7, height = 6.7 / 1.618 )

ggplot( model_power, aes( x = power, group = method, linetype = method ) ) + stat_ecdf( ) + facet_grid( . ~ heritability, scales = "free_x" )
dev.off( )
