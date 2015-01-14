library( ggplot2 )

argv = commandArgs( trailingOnly = TRUE )
if ( length( argv ) != 3 )
{
    message( "error: Wrong number of arguments." )
    message( "Usage: plot_roc title prediction output_file" )
    quit( )
}

title = argv[ 1 ]
prediction_file = argv[ 2 ]
output_file = argv[ 3 ]

pred = read.table( prediction_file, header = FALSE )
names( pred ) = c( "snp1", "snp2", "method", "score", "is_case" )
pred$tpr = rep( 0, length( pred$score ) )
pred$fpr = rep( 0, length( pred$score ) )

for( method in unique( pred$method ) )
{
    scores = pred$score[ pred$method == method ]
    correct = pred$is_case[ pred$method == method ]
    
    s_order = order( scores, decreasing = TRUE )
    pred$tpr[ pred$method == method ] = cumsum( correct[ s_order ] ) / sum( correct )
    pred$fpr[ pred$method == method ] = cumsum( 1 - correct[ s_order ] ) / sum( 1 - correct )
}

pdf( output_file, width = 8, height = 6 )

ggplot( pred ) + geom_line( aes( x = fpr, y = tpr, col = method ) ) +
    scale_x_continuous( "False positive rate") +
    scale_y_continuous( "True postive rate" ) +
    ggtitle( title ) +
    geom_abline( intercept = 0, slope = 1 )

dev.off( )
