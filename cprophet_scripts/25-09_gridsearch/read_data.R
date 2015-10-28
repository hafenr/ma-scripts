require(data.table)

readStatsData <- function(stats.file) {
    # Read output from gridsearch
    stats <- fread(stats.file)

    # Remove those parameter sets that scored no complex at all
    # so that no NaN will be produced when computing the decoy rate
    stats <- subset(stats, n_complexes_target_apmw_ok > 0)

    ## Compute rank based on number of identified target complexes
    # Compute 'apmw ok' version of decoy rate
    stats$complex_decoy_rate_apmw_ok <-
        stats$n_complexes_decoy_apmw_ok / stats$n_complexes_target_apmw_ok

    cat('Some parameter sets have a decoy rate of 0, these have to be removed. ',
        'To ensure that no parameter sets with a potentially very high score ',
        'are removed, there are their values for n_complexes_target_apmw_ok:\n',
        paste(stats[complex_decoy_rate_apmw_ok == 0,
              n_complexes_target_apmw_ok], collapse=' '))

    stats <- subset(stats, n_complexes_decoy_apmw_ok > 0)

    stats$diff_target_decoy_apwm_ok <-
        stats$n_complexes_target_apmw_ok - stats$n_complexes_decoy_apmw_ok
    stats$score_apmw <- 1 / stats$complex_decoy_rate_apmw_ok * stats$diff_target_decoy_apwm_ok
    stats$score <- 1 / stats$complex_decoy_rate * stats$diff_target_decoy

    ## Order the rows of the data set according to the target variable
    target.var.order <- order(stats$score_apmw, decreasing=TRUE)
    stats <- stats[target.var.order, ]

    stats
}

getCorumData <- function() {
    # Read corum identifiers
    corum <- fread('~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv', sep='\t',
                   stringsAsFactors=FALSE, colClasses=c(complex_id='character'))
    corum
}
