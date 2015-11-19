require(devtools)
require(ggplot2)
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

cutoffs <- seq(
    min(detected.features$apex_apmw_fit),
    max(detected.features$apex_apmw_fit),
    length=10
)

res <- makeROC(detected.features, manual.annotations, cutoffs)
plotROC(res)

# p <- p +
#     scale_x_continuous(limits=c(0, 1)) +
#     scale_y_continuous(limits=c(0, 1))
# print(p)

# abline(a=0, b=1, lty=2)

# tp.rate = sapply(cutoffs, calcTFR)
# fp.rate = FP / N = FP / (TN + FP)  

# calcRecall <- function(assessed.feats, )


# all.features <- assessComplexFeatures(
#     true.positive.features=manual.annotations,
#     detected.features=detected.features[apex_apmw_fit < 10000, ],
#     min.rt=3,
#     max.rt=84,
#     feature.vicinity.tol=5
# )
# n.points = all.features[type == 'FN', .N] + 
#            all.features[type == 'TP', .N] + 
#            all.features[type == 'FP', .N] + 
#            all.features[type == 'TN', .N] 
# stopifnot(n.points == nrow(all.features))

