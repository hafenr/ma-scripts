plotVennDiagram <- function(stats, n.topx=3, with.apmw=FALSE) {
    if (with.apmw) {
        ident.complexes <- strsplit(as.character(stats$target_complex_ids_apwm_ok), ',')
        targetvar.order <- order(stats$score_apmw, decreasing=T)
    } else {
        ident.complexes <- strsplit(as.character(stats$target_complex_ids), ',')
        targetvar.order <- order(stats$score, decreasing=T)
    }
    ## Plot overlap between detected corum complexes
    # Select the detected complex sets of the X parameter sets that have the highest
    # difference between target and decoys.
    stats.topx <- stats[targetvar.order[1:n.topx], ]
    ident.complexes.topx <- ident.complexes[targetvar.order[1:n.topx]]
    # Add the iteration number as the list's names (will later be used as the set name)
    names(ident.complexes.topx) <- as.character(targetvar.order[1:n.topx]-1)
    # Add all the corum complexes as a last element
    ident.complexes.topx$CORUM <- corum.complexes
    # Convert the list to a long format df
    dat <- melt(ident.complexes.topx, value.name='elements', variable.name='sets')
    v <- venneuler(dat)
    plot(v)
}


## Plot overlap between detected corum complexes
# Select the detected complex sets of the X parameter sets that have the highest
# difference between target and decoys.
plotVennDiagram(stats, n.topx=3, with.apmw=target.var == 'score_apmw')

