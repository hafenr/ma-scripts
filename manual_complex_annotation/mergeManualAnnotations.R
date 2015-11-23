require(devtools)
require(ggplot2)
require(hash)
# MACode git commit: f80708bd55d2d7724573a4d56599d078ce27e736
devtools::load_all('~/Dev/MACode')


annotations.1.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/151109_SEC-SWATH_CORUM_manual_annotation_mhe.tsv'
annotations.2.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/COMPLEXES_4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED_ANNOTATED.tsv'

manual.annotations.1.raw <- readManualAnnotationFile(annotations.1.raw)
manual.annotations.2.raw <- readManualAnnotationFile(annotations.2.raw)

complex.completeness.info <- manual.annotations.1.raw$complexes
complete.complex.ids <- complex.completeness.info[
    n_proteins_in_complex == n_proteins_in_complete_complex,
    complex_id
]

manual.annotations.1 <-
    createManualComplexAnnotations(manual.annotations.1.raw$annotations, 
                                   'apexes_fully_observed') 
manual.annotations.2 <-
    createManualComplexAnnotations(manual.annotations.2.raw$annotations,
                                   'apexes_fully_observed') 

manual.annotations <- mergeManualComplexAnnotations(manual.annotations.1,
                                                    manual.annotations.2)

manual.annotations.complete <-
    manual.annotations[complex_id %in% complete.complex.ids, ]

protein.traces.with.complex <- produceComplexAnnotatedProteinTraces(
    '~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv',
    '~/Dev/MAScripts/data/4_Fixed_weights_wMS1/Subsetted_OSW_Output/Peptide_Matrices/Long_Lists/4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED.tsv'
)
theoretically.possible.complexes <- unique(protein.traces.with.complex$complex_id)

# sec_complexes output form cprophet
detected.features <-
    fread('~/Dev/cprophet/run-674/iteration-0000/sec_complexes.tsv')

cutoffs <- seq(
    min(detected.features$apex_apmw_fit),
    max(detected.features$apex_apmw_fit) + 1,
    length=10
)

res <- makeROCWithSEC(detected.features, manual.annotations.complete, cutoffs)
plotROC(res, F)

# p <- p +
#     scale_x_continuous(limits=c(0, 1)) +
#     scale_y_continuous(limits=c(0, 1))
# print(p)

# abline(a=0, b=1, lty=2)

# tp.rate = sapply(cutoffs, calcTFR)
# fp.rate = FP / N = FP / (TN + FP)  

# calcRecall <- function(assessed.feats, )


# plotStatsForCutoff <- function(cval) {
#     filtered.feats = detected.features[apex_apmw_fit < cval, ]
#     all.features <- assessComplexFeatures(
#         true.positive.features=manual.annotations,
#         detected.features=filtered.feats,
#         min.rt=3,
#         max.rt=84,
#         feature.vicinity.tol=5
#     )
#     fn <- all.features[type == 'FN', .N]
#     tp <- all.features[type == 'TP', .N]
#     fp <- all.features[type == 'FP', .N]
#     tn <- all.features[type == 'TN', .N]
#     n.possible.features <- length(theoretically.possible.complexes) * length(seq(3, 84))
#     cat('CUTOFF: ', cval, '\n')
#     cat('#feats:\t', nrow(filtered.feats), '\n')
#     cat('SUM:\t', n.possible.features, '\n')
#     cat('SUM(is):\t', fn + tp + fp + tn, '\n')
#     cat('FNs:\t', fn, '\n')
#     cat('TPs:\t', tp, '\n')
#     cat('FPs:\t', fp, '\n')
#     cat('TNs:\t', tn, '\n')
# }
# plotStatsForCutoff(res$cutoff[9])
# plotStatsForCutoff(res$cutoff[10])


# Plot CUTOFF vs #TARGETS and CUTOFF vs #DECOYS
getNDecoys <- function(detected.features, cutoff) {
    feat.filtered <- detected.features[apex_apmw_fit < cutoff, ]
    feat.filtered[, is_decoy := grepl('^DECOY', complex_id)]
    length(unique(feat.filtered[is_decoy == T, complex_id]))
}
getNTargets <- function(detected.features, cutoff) {
    feat.filtered <- detected.features[apex_apmw_fit < cutoff, ]
    feat.filtered[, is_decoy := grepl('^DECOY', complex_id)]
    length(unique(feat.filtered[is_decoy == F, complex_id]))
}
n.decoys <- sapply(cutoffs, function(cutoff) getNDecoys(detected.features, cutoff))
n.targets <- sapply(cutoffs, function(cutoff) getNTargets(detected.features, cutoff))
decoy.target.data <-
    rbind(data.frame(count=n.decoys, type='decoy', cutoff=cutoffs),
          data.frame(count=n.targets, type='target', cutoff=cutoffs),
          data.frame(count=n.targets - n.decoys, type='target - decoy', cutoff=cutoffs))
p <- ggplot(decoy.target.data) +
    geom_line(aes(x=cutoff, y=count, color=type)) +
    geom_point(aes(x=cutoff, y=count, color=type))
print(p)


offsets <- seq(min(c(detected.features$left_apmw_fit,
                   detected.features$right_apmw_fit)),
               max(c(detected.features$left_apmw_fit,
                   detected.features$right_apmw_fit)),
               length=100)

assessComplexFeatures(manual.annotations,
                      detected.features,
                      2 * length(theoretically.possible.complexes))
roc <- makeROC(detected.features, manual.annotations.complete,
               cutoffs,
               2 * length(theoretically.possible.complexes))
plotROC(roc, F)

roc.data <- melt(roc, id.vars='cutoff', variable.name='type')
ggplot(roc.data) +
    geom_line(aes(x=cutoff, y=value, color=type))


# Plot FNs, TNs, FPs, TPs. as a function of cutoff in same plot

res <- sapply(cutoffs, function(cval) {
    assessComplexFeatures(manual.annotations,
                          detected.features[apex_apmw_fit < cval, ],
                          2 * length(theoretically.possible.complexes))
})
res <- as.data.frame(t(res))
res$cutoff <- cutoffs
res <- melt(res, id.var='cutoff', variable.name='type')
ggplot(aes(x=cutoff, y=value, color=type), data=res) +
    geom_line() + geom_point()


res <- sapply(offsets, function(offset) {
    assessComplexFeatures(
        manual.annotations,
        detected.features[
            left_apmw_fit >= 0 + offset & right_apmw_fit >= 0 + offset, 
        ],
        2 * length(theoretically.possible.complexes))
})
res <- as.data.frame(t(res))
res$offset <- offsets
res <- melt(res, id.var='offset', variable.name='type')
ggplot(aes(x=offset, y=value, color=type), data=res) +
    geom_line() + geom_point()










# SLIDING WINDOW CORRELATION CHECK
# CHECK HOW WELL THIS METHOD SEPARATES THE DECOY FROM REAL COMPLEXES

protein.traces.with.complex <- produceComplexAnnotatedProteinTraces(
    '~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv',
    '~/Dev/MAScripts/data/4_Fixed_weights_wMS1/Subsetted_OSW_Output/Peptide_Matrices/Long_Lists/4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED.tsv'
)
complex.protein.assoc <- unique(protein.traces.with.complex[, list(complex_id, protein_id)])
complex.protein.assoc[, complex_id := sample(complex_id)]


decoy.complex.assoc <- fread('~/Dev/MAScripts/data/decoy_complexes_seed_42.tsv')
intensity.info <- protein.traces.with.complex[, list(protein_id, sec, intensity)]
setkey(intensity.info)
intensity.info <- unique(intensity.info)
setkey(intensity.info, protein_id)
setkey(decoy.complex.assoc, protein_id)
protein.traces.decoy <- merge(decoy.complex.assoc, intensity.info,
                              allow.cartesian=T)

# protein.traces.decoy <- complex.protein.assoc[
#     subset(protein.traces.with.complex, select=-complex_id)
# ]
# protein.traces.decoy[, complex_id := paste0('DECOY_', complex_id)]
# protein.traces.decoy.subs <- portein[complex_id == 'DECOY_5342', ]
# ggplot(aes(x=sec, y=intensity, color=protein_id), data=protein.traces.decoy.subs) +
#     geom_line()
# ggplot(aes(x=sec, y=intensity, color=protein_id),
#        data=protein.traces[complex_id == '5615', ]) +
#     geom_line()

protein.traces <- rbind(
    protein.traces.with.complex[, list(complex_id, complex_name, protein_id, sec, intensity)],
    protein.traces.decoy[, list(complex_id, complex_name, protein_id, sec, intensity)]
)
window.size <- 15
# # TODO: protein.traces.with.complex should also have the decoy complexes in
# # them!
all.complex.ids <- unique(protein.traces$complex_id)
decoy.complex.ids <- unique(protein.traces.decoy$complex_id)
# target.complex.ids <- unique(protein.traces.with.complex$complex_id)
setkey(protein.traces, complex_id)

in.sample <- sapply(all.complex.ids, function(complex.id) {
    print(complex.id)
    protein.traces.subs <- protein.traces[complex_id == complex.id, ]
    # protein.traces.subs[, n_proteins_in_complete_complex := NULL]
    protein.traces.subs.wide <-
        dcast(protein.traces.subs, complex_id + protein_id + complex_name ~ sec,
              value.var='intensity')
    mat <- as.matrix(subset(protein.traces.subs.wide, select=-c(1, 2, 3)))
    if (nrow(mat) < 2) return(F)
    corr <- slidingWindowCorrelation(mat, window.size)
    print(corr)
    if (any(is.na(corr))) {
        F
    } else if (any(corr > 0.9)) {
        T
    } else {
        F
    }
    # corr.long <- data.table(corr=corr, sec=seq(min(protein.traces.subs$sec),
    #                                            max(protein.traces.subs$sec)),
    #                         complex_id=complex.id)
})

# We identify
# > sum(unique(manual.annotations$complex_id) %in% all.complex.ids[in.sample])
# [1] 425
# of
# > length(unique(manual.annotations$complex_id))
# [1] 591
# manually annotated complexes.
#
# For the decoys we have:
# > sum(grepl('^DECOY_', all.complex.ids[in.sample]))
# [1] 164
# of:
# > sum(grepl('^DECOY_', all.complex.ids))
# [1] 1888

max.intensity <- max(protein.traces.subs$intensity)
p <- ggplot(protein.traces.subs) +
            geom_point(aes(x=sec, y=intensity, color=protein_id)) +
            geom_line(aes(x=sec, y=intensity, color=protein_id)) +
            geom_line(aes(x=sec, y=corr * max.intensity/2), color='black',
                      data=corr.long) +
            geom_abline(intercept=max.intensity/2, linetype=2)

print(p)
