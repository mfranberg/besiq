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
names( model_power ) = c( "method", "heritability", "power", "lower", "upper", "maf", "sample_size" )

output_file = argv[ 5 ]

pdf( output_file, width = 2 * 6.7, height = 2 * 6.7 / 1.618 )

p = ggplot( model_power, aes( heritability, power, linetype = method ) ) + geom_line( ) + facet_grid( sample_size ~ maf, scales = "free_x" ) +
    scale_x_continuous( xlabel ) +
    scale_y_continuous( ylabel, limits = c( 0.0, 1.0 ) )

z <- ggplot_gtable(ggplot_build(p))

# add label for right strip
z <- gtable_add_cols(z, z$widths[[9]], 9)
z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, fill = gray(0.5))),
                          textGrob("Sample size", rot = -90, gp = gpar(col = gray(1)))),
                     4, 10, 8, name = paste(runif(2)))

# add label for top strip
z <- gtable_add_rows(z, z$heights[[3]], 2)
z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, fill = gray(0.5))),
                          textGrob("Minor allele frequency", gp = gpar(col = gray(1)))),
                     3, 4, 3, 8, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 9)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)

# draw it
# grid.newpage()
grid.draw(z)

dev.off( )
