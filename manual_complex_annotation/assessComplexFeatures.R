#' Check for each feature in `detected.features` if there is another feature in
#' `true.features` that is close to it. The endresult is a DF of feature
#' retention times that are flagged with either 'FP', 'TP', or 'FN'.
#' 
#' @param true.features data.table of true, manually annotated features.
#' The DF must have the columns: 'complex_id', 'rt'.
#' @param detected.features data.table of detected features.
#' The DF must have the columns: 'complex_id', 'rt'.
#' @param feature.vicinity.tol A number that indicates how close an
#' experimentally determined feature has to be to a manually annotated one, to
#' still count as a true positive.
#' @return A data.table with the columns: 'complex_id', 'rt', 'type', where
#' type is of type character an is either 'FP', 'FN', 'TP'.
assessComplexFeatures <- function(true.features, detected.features,
                                  feature.vicinity.tol=5) {
    complex.ids <- unique(detected.features$complex_id)

    do.call(rbind, lapply(complex.ids, function(cid) {
        rt.true <- true.features[complex_id == cid, ]$rt
        rt.exp <- detected.features[complex_id == cid, ]$center_rt
        # Declare integer vectors of true positives/false positives/false
        # negatives that are build up in the following loops.
        TPs <- integer(0)
        FPs <- integer(0)
        FNs <- integer(0)
        # Check for each experimentally found feature, i.e. RT...
        for (t.exp in rt.exp) {
            most.proximate.true.rt <- integer(0)
            smallest.delta.encountered <- Inf
            # Is the experimental rt close to some of the annotated ones?
            for(t.true in rt.true) {
                t.delta <- abs(t.true - t.exp)
                if (t.delta <= feature.vicinity.tol
                        && t.delta < smallest.delta.encountered) {
                    # This is the closest annotated RT to the experimental one.
                    # Save it.
                    most.proximate.true.rt <- t.true
                }
            }
            # Check if there was any annotated value that could be assigned to this
            # experimental rt.
            no.corresponding.true.rt.found <- length(most.proximate.true.rt) == 0
            if (no.corresponding.true.rt.found) {
                # No, this must be a false positive.
                FPs <- c(FPs, t.exp)
            } else {
                # Yes, this must be a true positive.
                TPs <- c(TPs, t.exp)
                # Remove the annotated vaue from the list so that it won't get
                # assigned to another feature rt of the same complex.
                rt.true <- setdiff(rt.true, most.proximate.true.rt)
            }
        }
        # All those annotated RT values that were not assigned (i.e. were also not
        # removed from the original array) by the setdiff call above are by
        # definition false negatives.
        FNs <- rt.true
        # The true negatives would be all theoretical complexes that aren't in TP.

        # Build a dataframe of the feature rts and add a character indicator
        # flag to what type they belong. 
        classifed.rts <- rbind(
            data.frame(rt=FNs, type=(if (length(FNs) > 0) 'FN' else character(0))),
            data.frame(rt=TPs, type=(if (length(TPs) > 0) 'TP' else character(0))),
            data.frame(rt=FPs, type=(if (length(FPs) > 0) 'FP' else character(0)))
        )

        classifed.rts$complex_id <- cid
        as.data.table(classifed.rts)
    }))
}


source('mergeManualAnnotations.R')
manual.annotations <- apex.dt.merged
# sec_complexes output form cprophet
detected.features <-
    fread('~/Dev/cprophet/run-674/iteration-0000/sec_complexes.tsv')
# Only look at the manual features were every measured 
# protein was also part of the peak (doesn't have to be
# a complete protein group).
true.features <- manual.annotations[apex_type == 'apexes_fully_observed']

assessComplexFeatures(true.features, detected.features, feature.vicinity.tol=5)
