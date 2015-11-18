require(devtools)
# MACode git commit: f80708bd55d2d7724573a4d56599d078ce27e736
devtools::load_all('~/Dev/MACode')


annotations.1.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/151109_SEC-SWATH_CORUM_manual_annotation_mhe.tsv'
annotations.2.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/COMPLEXES_4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED_ANNOTATED.tsv'

manual.annotations.1.raw <- readManualAnnotationFile(annotations.1.raw)
manual.annotations.2.raw <- readManualAnnotationFile(annotations.2.raw)

manual.annotations.1 <-
    createManualComplexAnnotations(manual.annotations.1.raw, 'apexes_fully_observed') 
manual.annotations.2 <-
    createManualComplexAnnotations(manual.annotations.2.raw, 'apexes_fully_observed') 

manual.annotations <- mergeManualComplexAnnotations(manual.annotations.1,
                                                    manual.annotations.2)

# sec_complexes output form cprophet
detected.features <-
    fread('~/Dev/cprophet/run-674/iteration-0000/sec_complexes.tsv')

# All features (detected as well as theoretically true) without the TNs
all.features <- assessComplexFeatures(
    true.positive.features=manual.annotations,
    detected.features=detected.features,
    feature.vicinity.tol=5
)

# Same as recall
calcTPR <- function(dt) {
    fn <- dt[type == 'FN', .N]
    tp <- dt[type == 'TP', .N]
    tp / (tp + fn)  # TP / P 
}
calcFPR <- function(dt) {
    fp <- dt[type == 'FP', .N]
    tn <- dt[type == 'TN', .N]
    fp / (tn + fp)  # FP / N
}

makeROC <- function(detected.features,
                    true.positive.features,
                    cutoffs,
                    feature.vicinity.tol=5) {
    fpr <- numeric(length=length(cutoffs))
    tpr <- numeric(length=length(cutoffs))
    for (i in seq_along(cutoffs)) {
        cat(sprintf('Calculating TPR/FPR. Iteration: %d\n', i))
        cval <- cutoffs[i]
        detected.features.filtered <-
            detected.features[apex_apmw_fit < cval, list(complex_id, center_rt)]
        if (nrow(detected.features.filtered) != 0) {
            all.features <- assessComplexFeatures(
                true.positive.features=true.positive.features,
                detected.features=detected.features.filtered,
                min.rt=3,
                max.rt=84,
                feature.vicinity.tol=5
            )
            fpr[i] <- calcFPR(all.features)
            tpr[i] <- calcTPR(all.features)
        } else {
            fpr[i] <- 0
            tpr[i] <- 0
        }
    }
    fpr.order <- order(fpr)
    list(FPR=fpr[fpr.order], TPR=tpr[fpr.order])
}

cutoffs <- seq(
    min(detected.features$apex_apmw_fit),
    max(detected.features$apex_apmw_fit),
    length=100
)
res <- makeROC(detected.features, manual.annotations, cutoffs)

plot(res$FPR[order(res$FPR)], res$TPR[order(res$FPR)],
     ylim=0:1,
     xlim=0:1,
     type='o',
     xlab='FPR | (1 - specificity)', ylab='TPR | sensitivity')
abline(a=0, b=1, lty=2)

# tp.rate = sapply(cutoffs, calcTFR)
# fp.rate = FP / N = FP / (TN + FP)  

# calcRecall <- function(assessed.feats, )



