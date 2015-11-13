require(data.table)
require(ggplot2)

readManualAnnotationFile <- function(fname) {
    annot <- fread(fname, sep='\t', sep2=',', colClasses=rep('character', 3))
    setnames(annot, c('complex_id', 'apexes_fully_observed',
                      'apexes_partially_observed'))
    annot
}

annotations.1.raw <-
    readManualAnnotationFile('151109_SEC-SWATH_CORUM_manual_annotation_mhe.tsv')
annotations.2.raw <- readManualAnnotationFile(
    'COMPLEXES_4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED_ANNOTATED.tsv'
)

# Split a list of apexes like '1,2,3,4' and create a data.frame
# data.frame(rt=c(1, 2, 3, 4), complex.id=X, apex.type=Y).
# If no split was possible, NULL is returned.
apexStringToDF <- function(complex.id, sep.apexes, apex.type) {
    rt <- as.numeric(strsplit(sep.apexes, ',')[[1]])
    if (length(rt) > 0) {  # if there was something to split on
        data.frame(complex_id=complex.id, rt=rt,
                   apex_type=apex.type,
                   stringsAsFactors=F)
    }
}

# Given a DF with columns 'complex_id', and another column holding 
# comma-separated strings of numbers, create a long list style DF
# where each row corresponds to a number.
createApexDF <- function(annotations, apex.col.name) {
    dframes.list <- mapply(apexStringToDF, annotations$complex_id,
                           annotations[[apex.col.name]], apex.col.name)
    # combine a list of dataframes into one large dataframe.
    # rbind will ignore entries that are NULL.
    apex.df <- do.call(rbind, dframes.list)
    rownames(apex.df) <- NULL
    apex.df
}
apex.df.1 <- rbind(createApexDF(annotations.1.raw, 'apexes_fully_observed'),
                   createApexDF(annotations.1.raw, 'apexes_partially_observed'))
apex.df.2 <- rbind(createApexDF(annotations.2.raw, 'apexes_fully_observed'),
                   createApexDF(annotations.2.raw, 'apexes_partially_observed'))

apex.dt.1 <- data.table(apex.df.1, key='complex_id')
apex.dt.2 <- data.table(apex.df.2, key='complex_id')

complex.ids <- unique(apex.dt.1$complex_id)

# Merge to vectors of numbers in such a way that the output won't
# contain numbers of the second vector that are within a interval
# [i - window, i + window] for each number i in the first vector.
# As an example, mergeRTs(c(1, 5), c(2, 3)) will result in c(1, 3, 5).
mergeRTs <- function(rts1, rts2, window=1) {
    ref.rts.with.spacings <- sapply(rts1, function(rt) {
        c(rt - window, rt, rt + window)
    })
    other.rts <- setdiff(rts2, ref.rts.with.spacings)
    merged.rts <- c(rts1, other.rts)
    merged.rts
}

stopifnot(mergeRTs(c(1, 5), c(3, 2)) == c(1, 5, 3))
stopifnot(mergeRTs(integer(0), c(3, 2)) == c(3, 2))
stopifnot(mergeRTs(c(3, 2), integer(0)) == c(3, 2))

# Merge the RTs for apexes of `apex.type` for two DTs.
# dt1 is treated as the reference DT.
createMergedList <- function(dt1, dt2, apex.type) {
    do.call(rbind, lapply(complex.ids, function(cid) {
        ref.rts <- dt1[complex_id == cid & apex_type == apex.type, rt]
        other.rts <- dt2[complex_id == cid & apex_type == apex.type, rt]
        merged <- mergeRTs(ref.rts, other.rts)
        if (length(merged) > 0) {
            data.table(complex_id=cid, rt=merged, apex_type=apex.type)
        } else {
            NULL
        }
    }))
}

apex.dt.merged <-
    rbind(createMergedList(apex.dt.1, apex.dt.2, 'apexes_fully_observed'),
          createMergedList(apex.dt.1, apex.dt.2, 'apexes_partially_observed'))
apex.dt.merged.otherway <-
    rbind(createMergedList(apex.dt.2, apex.dt.1, 'apexes_fully_observed'),
          createMergedList(apex.dt.2, apex.dt.1, 'apexes_partially_observed'))

ggplot(apex.dt.merged) +
    geom_density(aes(x=rt, fill=apex_type), alpha=0.5)
ggplot(apex.dt.merged.otherway) +
    geom_density(aes(x=rt, fill=apex_type), alpha=0.5)
