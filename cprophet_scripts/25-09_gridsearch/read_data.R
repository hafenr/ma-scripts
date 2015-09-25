readStatsData <- function(stats.file) {
    # Read output from gridsearch
    stats <- read.csv(stats.file, sep='\t')

    ## Compute rank based on number of identified target complexes
    # Compute 'apmw ok' version of decoy rate
    stats$complex_decoy_rate_apmw_ok <- stats$n_complexes_decoy_apmw_ok /
                                            stats$n_complexes_target_apmw_ok

    ## Impute missing/zero/infinite values since they cause problems downstream
    # Remove zero values from decoy rate
    is.zero <- stats$complex_decoy_rate_apmw_ok == 0
    # Impute Inf/NaN (source: 0/0 == NaN, x/0 == Inf)
    row.is.na <- is.na(stats$complex_decoy_rate_apmw_ok)
    is.inf <- is.infinite(stats$complex_decoy_rate_apmw_ok)
    stats <- subset(stats, !row.is.na & !is.inf & !is.zero)

    stats$diff_target_decoy_apwm_ok <- stats$n_complexes_target_apmw_ok -
                                        stats$n_complexes_decoy_apmw_ok
    stats$score_apmw <- 1 / stats$complex_decoy_rate_apmw_ok * stats$diff_target_decoy_apwm_ok
    stats$score <- 1 / stats$complex_decoy_rate * stats$diff_target_decoy

    ## Order the rows of the data set according to the target variable
    target.var.order <- order(stats$score_apmw, decreasing=TRUE)
    stats <- stats[target.var.order, ]

    stats
}

getCorumData <- function() {
    # Read corum identifiers
    corum <- read.csv('~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv', sep='\t',
                      stringsAsFactors=FALSE)
    corum
}
