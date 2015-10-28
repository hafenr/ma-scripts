require(data.table)
require(ggplot2)
require(prodlim)  # for row.match


# The file holding all peptide traces in long format.
# The file is downloadable from:
# 'https://drive.google.com/file/d/0B6_ub_yyGLS1bERZOXIzdGlQelk/view?usp=sharing'
peptr.file <- '~/Dev/MAScripts/data/peptide_traces_with_complex_id.tsv'

# Best performing: 1633
# Same but with remove_overlapping_features = False: 1291
run.name <- 'openms-complex-feat2'
iternum <- 674


## Get peptide traces
if (!exists('protein.traces')) {
    peptide.traces <- fread(peptr.file)
    # Sum peptide traces together to produce the protein traces
    protein.traces <- setnames(peptide.traces[, sum(intensity),
                               by=list(protein_id, complex_id, sec)],
                      'V1', 'intensity')
}

# decoy.complexes <- subset(fread('~/Dev/MAScripts/data/decoy_complexes_seed_42.tsv'),
#                           select=-complex_name)
# protein.traces.decoy <- subset(copy(protein.traces), select=-complex_id)
# setkey(decoy.complexes, protein_id)
# setkey(protein.traces.decoy, protein_id)
# protein.traces.annot.decoy <-
#     merge(decoy.complexes, protein.traces.decoy, allow.cartesian=T)

## Get the features that belong to proteins of identified complexes.
if (!exists('complex.feats.raw')) {
    complex.file <- file.path('/Volumes/rhafen/grid-search-runs',
                              run.name,
                              paste0('iteration-', sprintf('%04d', iternum)),
                              'sec_complexes.tsv')
    complex.feats.raw <- fread(complex.file)
    complex.feats.raw[, n_features:=length(center_rt), by=complex_id]

    # Add two new columns to differentiate between decoys and see which features
    # passed the apmw_ok test.
    complex.feats.raw$is_decoy <- grepl('^DECOY_', complex.feats.raw$complex_id)
    complex.feats.raw$apmw_ok <- complex.feats.raw$left_apmw_fit >= 0 & complex.feats.raw$right_apmw_fit >= 0
}

# Annotate each intensity observation with the complex it could theoretically
# be associated.
setkey(protein.traces, complex_id)
setkey(complex.feats.raw, complex_id)
# IMPORTANT: Decoy traces were removed in this step since they are not in the
# original peptide trace file
protein.traces.annot <- merge(protein.traces, complex.feats.raw, allow.cartesian=T)

# Produce different subset of features
complex.feats.target <- complex.feats.raw[is_decoy == F, ]
complex.feats.decoy <- complex.feats.raw[is_decoy == T, ]
complex.feats.target.apmwok <- complex.feats.target[is_decoy == F & apmw_ok == T, ]
complex.feats.decoy.apmwok <- complex.feats.raw[is_decoy == T & apmw_ok == T, ]


# ******************************************************************************

## Set the feature DT that will be used in the subsequent code
complex.feats <- complex.feats.target.apmwok

# ******************************************************************************


# Sample X complexes so the final plot won't be too large
detected.complex.ids <- unique(complex.feats$complex_id)
sample.complex.ids <- sample(detected.complex.ids, 20)
complex.feats.subset.pre <- complex.feats[complex_id %in% sample.complex.ids, ]
# compl <- compl[order(compl$d_score, decreasing=T), ]
# compl <- compl[compl$RT > 0, ]
# compl <- compl[n_features > 5, ]
# compl <- compl[1:10, ]
# compl <- complex.feats.subsetample(10), ]


# Compute the maximum intensity for the protein trace group (just for producing
# nicer plots).
max.int <- setnames(protein.traces.annot[, max(intensity), by=complex_id], 'V1', 'max_intensity')
setkey(max.int, complex_id)
setkey(complex.feats.subset.pre, complex_id)
complex.feats.subset <- merge(complex.feats.subset.pre, max.int)

# Subset the protein trace DT so it only contains observations that belong to
# complexes in the feature sample.

# Hack for deleting duplicate features
setkey(complex.feats.subset, NULL)
complex.feats.subset <- unique(complex.feats.subset)

protein.traces.subset <- protein.traces.annot[
    protein.traces.annot$complex_id %in% complex.feats.subset$complex_id, ]


p <- ggplot(protein.traces.subset)
p <- p + geom_line(aes(x=sec, y=intensity, color=protein_id))
p <- p + facet_wrap(~ complex_id, ncol=2, scales='free_y')
p <- p + geom_rect(aes(xmin=left_boundary_rt, xmax=right_boundary_rt,
                       ymin=0, ymax=max_intensity),
                   color='gray',
                   linetype='longdash',
                   alpha=0.1,
                   data=complex.feats.subset)
p <- p + geom_rect(aes(xmin=center_rt, xmax=center_rt, ymin=0, ymax=max_intensity),
                    # linetype='longdash',
                    color='#777777',
                    size=0.50,
                    data=complex.feats.subset)
# TODO: Add complex name to title
p <- p + theme(legend.position='none') +
         ggtitle(paste0('run: ', run.name, ' #', iternum))
plot(p)


# Find simular parameter sets
# params <- getParams(stats)
# params <- transform(
#     params,
#     param__remove_overlapping_peaks=as.character(param__remove_overlapping_peaks)
# )
# param.best <- params[1, ]
# query <- param.best
# query$param__remove_overlapping_peaks <- 'False'
# row.match(query, params)
# stats[107, 'iternum']

# d <- peptr.annot[complex_id == '157' | complex_id == '167', ]
# ggplot(d) +
#     geom_line(aes(x=sec, y=intensity, color=peptide_id)) +
#     facet_wrap(~ complex_id) +
#     theme(legend.position='none')


# ggplot(complex.feats.raw) +
#     geom_density(aes(x=RT, color=apmw_ok))
