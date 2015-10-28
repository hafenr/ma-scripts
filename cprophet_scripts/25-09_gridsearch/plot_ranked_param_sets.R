source('read_data.R')

## Plot number of targets and decoys for parameter sets ranked by their target score
plotRankedParameterSets <- function(stats, target.var, xlim=c(0, nrow(stats))) {
    target.var.order <- order(stats[[target.var]], decreasing=TRUE)
    plot(target.var.order, stats[[target.var]], xlim=xlim)
    points(target.var.order, stats$n_complexes_decoy_apmw_ok, col='red')
    points(target.var.order, stats$n_complexes_target_apmw_ok, col='green', pch=3)
}

plotRankedParameterSets(stats, 'score_apmw')
