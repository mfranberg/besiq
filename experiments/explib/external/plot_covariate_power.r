library( ggplot2 )

argv = commandArgs( trailingOnly = TRUE )
if ( length( argv ) != 5 )
{
    message( "error: Wrong number of arguments." )
    message( "Usage: plot_model_power xlabel ylabel title cov_model output_file" )
    quit( )
}

model_power = read.table( argv[ 4 ], header = F )
names( model_power ) = c( "method", "effect", "power", "lower", "upper", "cov" )
model_power$cov[ model_power$cov == 0 ] = "Without covariates"
model_power$cov[ model_power$cov == 1 ] = "With covariates"

xlabel = argv[ 1 ]
ylabel = argv[ 2 ]
title = argv[ 3 ]
output_file = argv[ 5 ]

pdf( output_file )

dodge = position_dodge( width = 0.4 )
ggplot( model_power, aes( x = effect, y = power, group = method, ymin = lower, ymax = upper ) ) +
    geom_errorbar( width = 0.2, position = dodge ) +
    geom_point( aes( shape = method, fill = cov ), position = dodge, size = 2.5 ) +
    scale_x_continuous( xlabel ) +
    scale_y_continuous( ylabel, limits = c( 0.0, 1.0 ) ) +
    ggtitle( title ) +
    scale_shape_manual( values = c( 21, 22, 23 ) ) +
    scale_fill_manual( values = c( "black", "white"  ) ) +
    guides( fill = guide_legend( override.aes = list( shape = 21 ) ) ) +
    theme_bw( )

dev.off( )
