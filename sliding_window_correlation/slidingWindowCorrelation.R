require(data.table)
require(ggplot2)

rm(list=ls())

complex.id <- '564-6'

produceComplexAnnotatedProteinTraces <- function(filename.corum.complex.assoc,
                                                 filename.peptide.traces.long) {
    # Read corum identifiers
    corum.protein.assoc <-
        fread(filename.corum.complex.assoc, sep='\t',
              stringsAsFactors=FALSE, colClasses=c(complex_id='character'))

    # Read peptide traces and produce protein trace
    peptide.traces <- fread(filename.peptide.traces.long)
    # Sum peptide traces together to produce the protein traces
    protein.traces <- setnames(peptide.traces[, sum(peptide_intensity),
                               by=list(protein_id, sec)],
                      'V1', 'intensity')

    # Add a new column to the corum proteins that tells how many proteins
    # are theoretically in the corum complex.
    corum.protein.assoc[
        ,
        n_proteins_in_complete_complex := length(unique(protein_id)),
        by=complex_id
    ]

    # Merge protein traces together with a list of corum complexes
    # so that each intensity observation is annotated with the complex id
    # it __could potentially__ belong to.
    setkey(protein.traces, protein_id)
    setkey(corum.protein.assoc, protein_id)
    protein.traces.with.complex <-
        merge(protein.traces, corum.protein.assoc, allow.cartesian=TRUE)

    protein.traces.with.complex
}

protein.traces.with.complex <- produceComplexAnnotatedProteinTraces(
    '~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv',
    '~/Dev/MAScripts/data/4_Fixed_weights_wMS1/Subsetted_OSW_Output/Peptide_Matrices/Long_Lists/4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED.tsv'
)

protein.traces.subs <- protein.traces.with.complex[complex_id == complex.id, ]
protein.traces.subs[, n_proteins_in_complete_complex := NULL]

protein.traces.subs.wide <-
    dcast(protein.traces.subs, complex_id + protein_id + complex_name ~ sec,
          value.var='intensity')
mat <- as.matrix(subset(protein.traces.subs.wide, select=-c(1, 2, 3)))

computeWindowedCorrelation <- function(trace.mat, start.window.idx, window.size) {
    # Make sure that the window matrix will actually be a matrix and not just
    # a 1-column vector.
    if (start.window.idx == ncol(trace.mat)) {
        start.window.idx <- ncol(trace.mat) - 1
    }
    end.window.idx <- min(ncol(trace.mat), start.window.idx + window.size)
    window.trace.mat <- trace.mat[, start.window.idx:end.window.idx]

    # Rows with only 0 will lead to uncomputable standard deviations when
    # computing the correlation. Those rows are imputed with Norm(0, 1) noise.
    is.all.zero.row <- rowSums(window.trace.mat) == 0
    if (any(is.all.zero.row)) {
        # window.trace.mat <- window.trace.mat[!is.all.zero.row, ]
        effective.window.size <- end.window.idx - start.window.idx + 1
        n.zero.entries <- effective.window.size * sum(is.all.zero.row)
        window.trace.mat[is.all.zero.row, ] <- rnorm(n.zero.entries)
    }

    # Compute the mean correlation between all the traces within this window.
    # The mean is computed only on the upper triangular matrix.
    correls <- cor(t(window.trace.mat), t(window.trace.mat))
    correls.upper <- correls[upper.tri(correls)]
    r <- mean(correls.upper)
    if (!is.na(r)) r else 0
}

window.size <- 10

corr <- sapply(seq(1, ncol(mat)), function(i) {
    computeWindowedCorrelation(mat, i, window.size)
})

corr.long <- data.table(corr=corr, sec=seq(min(protein.traces.subs$sec),
                                           max(protein.traces.subs$sec)),
                        complex_id=complex.id)

max.intensity <- max(protein.traces.subs$intensity)
p <- ggplot(protein.traces.subs) +
            geom_point(aes(x=sec, y=intensity, color=protein_id)) +
            geom_line(aes(x=sec, y=intensity, color=protein_id)) +
            geom_line(aes(x=sec, y=corr * max.intensity/2), color='black',
                      data=corr.long) +
            geom_abline(intercept=max.intensity/2, linetype=2)

print(p)


