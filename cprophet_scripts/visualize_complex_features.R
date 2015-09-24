require(data.table)
require(ggplot2)
require(prodlim)  # for row.match

# rm(list=ls())

if (!exists('peptr')) {
    # Downloadable from:
    # 'https://drive.google.com/file/d/0B6_ub_yyGLS1bERZOXIzdGlQelk/view?usp=sharing'
    file.loc <- '~/Dev/MAScripts/data/peptide_traces_with_complex_id.tsv'
    peptr <- fread(file.loc)
    peptr <- setnames(peptr[, sum(intensity), by=list(protein_id, sec)],
                      'V1', 'intensity')
}

# Best performing: 1633
# Same but with remove_overlapping_features = False: 1291
run.name <- 'openms'
iternum <- 1226
if (!exists('compl.raw')) {
  complex.file <- file.path('/Volumes/rhafen/grid-search-runs',
                            run.name,
                            paste0('iteration-', iternum),
                            'sec_complexes.tsv')
  compl.raw <- fread(complex.file)
  compl.raw[, n_features:=length(RT), by=complex_id]

  is.decoy <- grepl('^DECOY_', compl.raw$complex_id)
  compl.raw$is_decoy <- is.decoy
  compl.raw$apmw_ok <- compl.raw$left_apmw_fit >= 0 & compl.raw$right_apmw_fit >= 0
}
compl.target <- compl.raw[!is.decoy, ]
compl.decoy <- compl.raw[is.decoy, ]
compl.decoy.apmwok <- compl.decoy[compl.decoy$apmw_ok, ]
compl.target.apmwok <- compl.target[compl.target$apmw_ok ]


compl <- compl.target.apmwok

detected.complex.ids <- unique(compl$complex_id)
sample.complex.ids <- sample(detected.complex.ids, 20)
compl <- compl[complex_id %in% sample.complex.ids, ]
# compl <- compl[order(compl$d_score, decreasing=T), ]
# compl <- compl[compl$RT > 0, ]
# compl <- compl[n_features > 5, ]
# compl <- compl[1:10, ]
# compl <- compl[sample(10), ]


setkey(peptr, protein_id)
setkey(compl, protein_id)

peptr.annot <- merge(peptr, compl)

max.int <- setnames(peptr.annot[, max(intensity), by=complex_id], 'V1', 'max_intensity')

setkey(max.int, complex_id)
setkey(compl, complex_id)
compl <- compl[max.int]

p <- ggplot(peptr.annot)
p <- p + geom_line(aes(x=sec, y=intensity, color=protein_id)) +
    facet_wrap(~ complex_id, ncol=2, scales='free_y') +
    geom_rect(aes(xmin=leftWidth, xmax=rightWidth,
                  ymin=0, ymax=max_intensity),
              # linetype='longdash',
              alpha=0.1,
              color='gray',
              data=compl)
p <- p + geom_rect(aes(xmin=RT, xmax=RT, ymin=0, ymax=max_intensity),
                    # linetype='longdash',
                    color='#777777',
                    size=0.50,
                    data=compl)
p <- p + theme(legend.position='none') +
         ggtitle(paste0('run: ', run.name, ' #', iternum))
plot(p)

# Find simular parameter sets
params <- getParams(stats)
params <- transform(
    params,
    param__remove_overlapping_peaks=as.character(param__remove_overlapping_peaks)
)
param.best <- params[1, ]
query <- param.best
query$param__remove_overlapping_peaks <- 'False'
row.match(query, params)
stats[107, 'iternum']

# d <- peptr.annot[complex_id == '157' | complex_id == '167', ]
# ggplot(d) +
#     geom_line(aes(x=sec, y=intensity, color=peptide_id)) +
#     facet_wrap(~ complex_id) +
#     theme(legend.position='none')


# ggplot(compl.raw) +
#     geom_density(aes(x=RT, color=apmw_ok))
