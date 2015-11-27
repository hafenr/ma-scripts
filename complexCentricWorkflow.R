require(devtools)
require(ggplot2)
devtools::load_all('~/Dev/MACode')


protein.traces.long <-
    produceComplexAnnotatedProteinTraces(e4.peptide.traces.wide.filtered,
                                         corum.complex.protein.assoc)

protein.traces.wide <- longProtTracesToWide(protein.traces.long)
protein.traces.wide <- protein.traces.wide[, complex_name := NULL]

setnames(protein.traces.wide, 'protein_id', 'peptide_id')
setnames(protein.traces.wide, 'complex_id', 'protein_id')

write.table(protein.traces.wide, '~/Desktop/protein_traces_wide_CC.tsv',
            sep='\t', row.names=F)

## now run cprophet like this:
# import cprophet
#
# cprophet.pipeline.run_pipeline(
#     '/Users/robin/Desktop/protein_traces_wide_CC.tsv',
#     '/Users/robin/Desktop/cprophet/input_proteins.csv',
#     '/Users/robin/Desktop/cprophet/input_sec.txt',
#     '/Users/robin/Desktop/cc_out/',
#     params=cprophet.pipeline.DEFAULT_PARAMS
# )

# Read detected features (these are now complex features)
complex.features <- fread('~/Desktop/cc_out_default_params/sec_proteins_with_dscore.csv')

# Convert to the same style of complex feature data.table that is used by the
# PC workflow. This enables us to use all other functions that operate on
# complex features.
detected.features <- convertToPCComplexFeatureFormat(complex.features) 

# Copy all score columns back into the complex features data.frame
score.column.names <- colnames(complex.features)[
    grep('^var_', colnames(complex.features))
]
for (col in c(score.column.names, 'd_score')) {
    detected.features[[col]] <- complex.features[[col]]
}

detected.features[, is_decoy := grepl('^DECOY', complex_id)]

annot <- annotateComplexFeatures(detected.features, manual.annotations.full)
# Add column is_true_positive back into dt
detected.features$is_true_positive <- annot$is_true_positive

plot(density(detected.features$d_score), main='d_score distribution',
     xlab='d_score', ylab='density')
lines(density(detected.features[is_true_positive == TRUE, d_score]), col='green')
lines(density(detected.features[is_decoy == TRUE, d_score]), col='red')
lines(density(detected.features[is_decoy == FALSE, d_score]), col='blue')
legend('topright', c('all', 'true positive', 'decoy', 'target'),
       col=c('black', 'green', 'red', 'blue'),
       lty=1, bg='#eeeeee')

# detected.features <- rbind(
#     detected.features,
#     detected.features[is_true_positive == T],
#     detected.features[is_true_positive == T],
#     detected.features[is_true_positive == T],
#     detected.features[is_true_positive == T]
# )

require(MASS)
lda.fit <- lda(is_true_positive ~ var_xcorr_coelution + var_intensity_score +
                                  var_log_sn_score + var_elution_model_fit_score,
               data=detected.features)
# cv.fit <- cvFit(
#     lda, formula=is_true_positive ~ var_xcorr_coelution + var_intensity_score +
#     var_log_sn_score + var_elution_model_fit_score,
#     data=detected.features,
#     cost=function(true, pred) { mean(true != pred) },
#     K=10
# )
# print(cv.fit)

theoretically.possible.complexes <- unique(protein.traces.long$complex_id)
params <- seq(min(complex.features$d_score), max(complex.features$d_score),
              length=100) 
roc <- makeROC(
        detected.features,   # detected feats
        manual.annotations,  # true pos
        params,
        # the number of complexes that were present in the input data
        2 * length(theoretically.possible.complexes),  
        function(feats, p) {
            feats[d_score > p, ]
        })
plotROC(roc, F)

protein.traces.wo.complex <- subset(protein.traces.long, select=-c(complex_name, complex_id))
setkey(protein.traces.wo.complex)
protein.traces.wo.complex <- unique(protein.traces.wo.complex)

pdf('~/Desktop/features_w_params_run_647.pdf')
complex.features.fil <- 
    complex.features.fil[order(complex.features.fil$d_score, decreasing=T), ]
for (i in 1:nrow(complex.features.fil)) {
    feat <- complex.features.fil[i, ]
    cat('DECOY: ', feat$decoy == 1, '\n')
    cat('D_SCORE: ', feat$d_score, '\n')
    subunits <- strsplit(feat$aggr_Fragment_Annotation, ';')[[1]]
    trace.data <- protein.traces.wo.complex[protein_id %in% subunits, ]
    complex.id <- gsub('pep_', '', feat$transition_group_id)
    feat$max_intensity <- max(trace.data$intensity)
    p <- ggplot(aes(x=sec, y=intensity, color=protein_id),
           data=trace.data) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept=feat$leftWidth) +
        geom_vline(xintercept=feat$rightWidth) +
        ggtitle(complex.id)
    max.sec <- max(trace.data$sec)
    max.intensity <- max(trace.data$intensity)
    p <- p + annotate('text', x=0.88*max.sec, y=0.96*max.intensity,
                      label=round(feat$d_score, 2), size=10)
    print(p)
    # readline(prompt = "Pause. Press <Enter> to continue...")
}
dev.off()


complex.ids.true <- unique(manual.annotations$complex_id)

# Plot CUTOFF vs #TARGETS and CUTOFF vs #DECOYS
nDecoys <- function(dt) {
    dt[, is_decoy := grepl('^DECOY', complex_id)]
    length(unique(dt[is_decoy == T, complex_id]))
}
nTrueTargets <- function(dt) {
    sum(complex.ids.true %in% dt[, complex_id])
}
nTargets <- function(dt) {
    dt[, is_decoy := grepl('^DECOY', complex_id)]
    length(unique(dt[is_decoy == F, complex_id]))
}

n.decoys <- sapply(params, function(p) {
    nDecoys(detected.features[ d_score > p, ])
})
n.targets <- sapply(params, function(p) {
    nTargets(detected.features[ d_score > p, ])
})
n.true.targets <- sapply(params, function(p) {
    nTrueTargets(detected.features[ d_score > p, ])
})

decoy.target.data <-
    rbind(data.frame(count=n.decoys, type='decoy', param=params),
          data.frame(count=n.targets, type='target', param=params),
          data.frame(count=n.true.targets, type='true targets', param=params)
          # , data.frame(count=n.targets - n.decoys, type='target - decoy',
          #            param=params)
          )
p <- ggplot(decoy.target.data) +
    geom_line(aes(x=param, y=count, color=type)) +
    geom_point(aes(x=param, y=count, color=type))
print(p)


