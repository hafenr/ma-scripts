require(devtools)
devtools::load_all('~/Dev/MACode')



# Compute protein traces and annotate them with a complex id.
# Format will be long list.
protein.traces.with.complex <- produceComplexAnnotatedProteinTraces(
    '~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv',
    '~/Dev/MAScripts/data/4_Fixed_weights_wMS1/Subsetted_OSW_Output/Peptide_Matrices/Long_Lists/4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED.tsv'
)

# Extract only those intensities that belong to some complex id(s).
complex.id <- '192'
protein.traces.subs <- protein.traces.with.complex[complex_id == complex.id, ]
protein.traces.subs[, n_proteins_in_complete_complex := NULL]

# Convert the long list to wide format and then to matrix.
protein.traces.subs.wide <-
    dcast(protein.traces.subs, complex_id + protein_id + complex_name ~ sec,
          value.var='intensity')
mat <- as.matrix(subset(protein.traces.subs.wide, select=-c(1, 2, 3)))

# Compute the sliding window correlation.
corr <- slidingWindowCorrelation(mat, window.size=15, score='diff')
# Convert back to long list.
corr.long <- data.table(corr=corr, sec=seq(min(protein.traces.subs$sec),
                                           max(protein.traces.subs$sec)),
                        complex_id=complex.id)
max.intensity <- max(protein.traces.subs$intensity)
p <- ggplot(protein.traces.subs) +
            geom_point(aes(x=sec, y=intensity, color=protein_id)) +
            geom_line(aes(x=sec, y=intensity, color=protein_id)) +
            geom_line(aes(x=sec, y=corr * max.intensity/2), color='black',
                      data=corr.long) +
            geom_abline(intercept=max.intensity/2, linetype=3) +
            geom_abline(intercept=max.intensity/2 * 0.9, linetype=3,
                        color='red') +
            ggtitle(paste0('CORUM complex ', complex.id))
print(p)
