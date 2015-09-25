require(data.table)
require(ggplot2)
require(prodlim)  # for row.match


# The file holding all peptide traces in long format.
# The file is downloadable from:
# 'https://drive.google.com/file/d/0B6_ub_yyGLS1bERZOXIzdGlQelk/view?usp=sharing'
peptr.file <- '~/Dev/MAScripts/data/peptide_traces_with_complex_id.tsv'

# Best performing: 1633
# Same but with remove_overlapping_features = False: 1291
run.name <- 'openms'
iternum <- 1633


## Get peptide traces
if (!exists('peptr')) {
    peptr <- fread(peptr.file)
    # Sum peptide traces together to produce the protein traces
    peptr <- setnames(peptr[, sum(intensity), by=list(protein_id, sec)],
                      'V1', 'intensity')
}

## Get the features that belong to proteins of identified complexes.
if (!exists('compl.raw')) {
  complex.file <- file.path('/Volumes/rhafen/grid-search-runs',
                            run.name,
                            paste0('iteration-', iternum),
                            'sec_complexes.tsv')
  compl.raw <- fread(complex.file)
  compl.raw[, n_features:=length(RT), by=complex_id]

  # Add two new columns to differentiate between decoys and see which features
  # passed the apmw_ok test.
  compl.raw$is_decoy <- grepl('^DECOY_', compl.raw$complex_id)
  compl.raw$apmw_ok <- compl.raw$left_apmw_fit >= 0 & compl.raw$right_apmw_fit >= 0
}

# Annotate each intensity observation with the complex it could theoretically
# be associated.
setkey(peptr, protein_id)
setkey(compl.raw, protein_id)
peptr.annot <- merge(peptr, compl.raw, allow.cartesian=T)

# Produce different subset of features
compl.target <- compl.raw[is_decoy == F, ]
compl.decoy <- compl.raw[is_decoy == T, ]
compl.target.apmwok <- compl.target[is_decoy == F & apmw_ok == T, ]
compl.decoy.apmwok <- compl.raw[is_decoy == T & apmw_ok == T, ]


# ******************************************************************************

## Set the feature DT that will be used in the subsequent code
compl <- compl.target.apmwok

# ******************************************************************************


# Sample X complexes so the final plot won't be too large
detected.complex.ids <- unique(compl$complex_id)
sample.complex.ids <- sample(detected.complex.ids, 20)
compl.s <- compl[complex_id %in% sample.complex.ids, ]
# compl <- compl[order(compl$d_score, decreasing=T), ]
# compl <- compl[compl$RT > 0, ]
# compl <- compl[n_features > 5, ]
# compl <- compl[1:10, ]
# compl <- compl[sample(10), ]


# Compute the maximum intensity for the protein trace group (just for producing
# nicer plots).
max.int <- setnames(peptr.annot[, max(intensity), by=complex_id], 'V1', 'max_intensity')
setkey(max.int, complex_id)
setkey(compl.s, complex_id)
compl.s <- merge(compl.s, max.int)

# Subset the protein trace DT so it only contains observations that belong to
# complexes in the feature sample.
peptr.s <- peptr.annot[peptr.annot$complex_id %in% compl.s$complex_id,  ]
p <- ggplot(peptr.s)
p <- p + geom_line(aes(x=sec, y=intensity, color=protein_id)) +
    facet_wrap(~ complex_id, ncol=2, scales='free_y') +
    geom_rect(aes(xmin=leftWidth, xmax=rightWidth,
                  ymin=0, ymax=max_intensity,
                  color=protein_id),
              linetype='longdash',
              alpha=0.1,
              data=compl.s)
p <- p + geom_rect(aes(xmin=RT, xmax=RT, ymin=0, ymax=max_intensity),
                    # linetype='longdash',
                    color='#777777',
                    size=0.50,
                    data=compl.s)
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


# ggplot(compl.raw) +
#     geom_density(aes(x=RT, color=apmw_ok))
