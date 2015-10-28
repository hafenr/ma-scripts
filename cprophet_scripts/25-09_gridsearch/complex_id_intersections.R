source('read_data.R')

corum.data <- getCorumData()
corum.complexes <- corum.data$complex_id
with.apmw <- TRUE
n.topx <- 4

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

# Compute how many complexes overlap between the different runs
n.topx <- 200
cat('Overlap between different param sets in #complexes:\n')
cat('(Numbers indicate position in top-scorer list)\n')
overlap.matrix <- matrix(0, nrow=n.topx, ncol=n.topx)
for (i in 1:n.topx) {
    for (j in 1:n.topx) {
        overlap <- length(intersect(ident.complexes[[i]],
                                    ident.complexes[[j]]))
        cat('|', i, '∩', j, '|', '=', overlap, '\n')
        overlap.matrix[i, j] <- overlap
    }
}
# overlap.matrix[!upper.tri(overlap.matrix)] <- 0
print(overlap.matrix)
cat('Overlap with CORUM:', '\n')
for (set.idx in 1:n.topx) {
    overlap <- length(intersect(ident.complexes[[set.idx]],
                                corum.complexes))
    cat('|', set.idx, '∩ CORUM', '|', '=', overlap, '\n')
}


# Identify core complexes
n.complexes <- 1:200
intersection.size <- sapply(n.complexes, function(i) {
    length(Reduce(intersect, ident.complexes[1:i]))
})
plot(n.complexes, intersection.size, type='l')



##############################################################################
#              How many times each complex occursin a paramset               #
##############################################################################
all.ident.complexes <- Reduce(union, ident.complexes)
complex.occurrence.counts <- sapply(all.ident.complexes, function(complex.id) {
    is.in.complex.sets <- sapply(ident.complexes, function(complex.set) {
        complex.id %in% complex.set
    })
    sum(is.in.complex.sets)
})
complex.occurrence.counts.df <- data.frame(
    complex_id=names(complex.occurrence.counts),
    count=complex.occurrence.counts
)

corum.complex.protein.assoc = getCorumData()
corum.complexes = unique(corum.complex.protein.assoc[, .(complex_id, complex_name)])
setkey(corum, complex_id)

complex.occurrence.counts.df <- merge(complex.occurrence.counts.df, corum.complexes)

ord <- order(complex.occurrence.counts.df$count, decreasing=T)
complex.occurrence.counts.df <- complex.occurrence.counts.df[ord, c(1, 3, 2)]
rownames(complex.occurrence.counts.df) <- 1:nrow(complex.occurrence.counts.df)
print(complex.occurrence.counts.df[1:20, ])


##############################################################################
#                     Identify core subset of complexes                      #
##############################################################################
# complex.core.set <- as.character(complex.occurrence.counts.df[1:50, 'complex_id'])
complex.core.set <- Reduce(intersect, ident.complexes[1:10])
complex.train.set <- Reduce(union, ident.complexes[1:100])

complex.core.set <- complex.train.set

peptr.file <- '~/Dev/MAScripts/data/peptide_traces_with_complex_id.tsv'
if (!exists('protein.traces')) {
    peptide.traces <- fread(peptr.file)
    # Sum peptide traces together to produce the protein traces
    protein.traces <- setnames(peptide.traces[, sum(intensity),
                               by=list(protein_id, complex_id, sec)],
                      'V1', 'intensity')
}
# Add name column
setkey(protein.traces, complex_id)
protein.traces.w.name <- protein.traces[corum.complexes]
protein.traces.core <- protein.traces.w.name[complex_id %in% complex.core.set, ]
protein.traces.core$extended_complex_id <-
    paste(protein.traces.core$complex_id,
          protein.traces.core$complex_name,
          sep=' - ')


# Merge in features of top run
run.name <- 'openms-complex-feat3'
iternum <- 674
complex.file <- file.path('/Volumes/rhafen/grid-search-runs',
                          run.name,
                          paste0('iteration-', sprintf('%04d', iternum)),
                          'sec_complexes.tsv')
complex.feats.raw <- fread(complex.file)
setkey(complex.feats.raw, complex_id)
setkey(protein.traces.core, complex_id)
protein.traces.core.w.feats <- merge(protein.traces.core, complex.feats.raw, by='complex_id', allow.cartesian=T)

# Produce a plot of them
pdf('union_top_100_runs_complexes.pdf')
for (complex.id in complex.core.set) {
    print(complex.id)
    traces <- protein.traces.core[complex_id == complex.id, ]
    features <- complex.feats.raw[complex_id == complex.id, ]
    complex.name <- traces[1, complex_name]
    p = ggplot(traces) +
        geom_line(aes(x=sec, y=intensity, color=protein_id)) +
        theme(legend.position='none') +
        ggtitle(sprintf('%s | %s', complex.id, complex.name))
    # p = p + geom_rect(aes(xmin=left_boundary_rt, xmax=right_boundary_rt,
    #                    ymin=0, ymax=approx_peak_height),
    #                color='gray',
    #                linetype='longdash',
    #                alpha=0.1,
    #               data=features) +
    #     geom_rect(aes(xmin=center_rt, xmax=center_rt, ymin=0,
    #                   ymax=approx_peak_height),
    #                 # linetype='longdash',
    #                 color='#777777',
    #                 size=0.50,
    #               data=features)
    print(p)
}
dev.off()
