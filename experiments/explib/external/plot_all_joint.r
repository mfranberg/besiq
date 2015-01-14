library( ggplot2 )
library( gtable )
library( proto )

### GGPLOT reversed ecdf

stat_ecdf_reversed <- function (mapping = NULL, data = NULL, geom = "step", position = "identity", n = NULL, ...) {
    StatEcdfReversed$new(mapping = mapping, data = data, geom = geom, position = position, n = n, ...)
}

StatEcdfReversed <- proto(ggplot2:::Stat, {
    objname <- "ecdf"
    
    calculate <- function(., data, scales, n = NULL, ...) {
        
        # If n is NULL, use raw values; otherwise interpolate
        if (is.null(n)) {
            xvals <- unique(data$x)
        } else {
            xvals <- seq(min(data$x), max(data$x), length.out = n)
        }
        
        y <- 1 - ecdf(data$x)(xvals)
        
        # make point with y = 0, from plot.stepfun
        rx <- range(xvals)
        if (length(xvals) > 1L) {
            dr <- max(0.08 * diff(rx), median(diff(xvals)))
        } else {
            dr <- abs(xvals)/16
        }
        
        x0 <- rx[1] - dr
        x1 <- rx[2] + dr
        y0 <- 1
        y1 <- 0
        
        data.frame(x = c(x0, xvals, x1), y = c(y0, y, y1))
    }
    
    default_aes <- function(.) aes(y = ..y..)
    required_aes <- c("x")
    default_geom <- function(.) ggplot2:::GeomStep
})

### End GGPLOT reversed ecdf

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
names( model_power ) = c( "Method", "model", "heritability", "power", "lower", "upper", "maf", "sample_size" )
# model_power = subset( model_power, heritability == 0.02 )

model_power$Method = factor( model_power$Method, levels( model_power$Method )[ c( 3, 2, 1, 4 ) ] )

output_file = argv[ 5 ]

pdf( output_file, width = 2 * 6.7, height = 2 * 6.7 / 1.618 )

#p = ggplot( model_power, aes( x = power, linetype = Method ) ) + stat_ecdf_reversed( geom = "smooth", colour = "black" ) + facet_grid( sample_size ~ maf, scales = "free_x" ) +
p = ggplot( model_power, aes( x = power, colour = Method ) ) + stat_ecdf_reversed( geom = "smooth" ) + facet_grid( sample_size ~ heritability, scales = "free_x" ) +
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
                          textGrob("Heritability", gp = gpar(col = gray(1)))),
                     3, 4, 3, 8, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 9)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)

# draw it
# grid.newpage()
grid.draw(z)

dev.off( )
