require('data.table')
require('ggplot2')

rm(list=ls())

# Read corum identifiers
corum.protein.assoc <-
    fread('~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv', sep='\t',
          stringsAsFactors=FALSE, colClasses=c(complex_id='character'))

# Read peptide traces and produce protein trace
peptide.traces.file <-
    '~/Dev/MAScripts/data/4_Fixed_weights_wMS1/Subsetted_OSW_Output/Peptide_Matrices/Long_Lists/4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED.tsv'
peptide.traces <- fread(peptide.traces.file)
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

# Add a new column to the protein traces which specifies how many
# proteins were actually merged to each corum complex.
# are theoretically in the corum complex.
protein.traces.with.complex[
    , n_proteins_in_complex := length(unique(protein_id)),
    by=complex_id
]

complex.ids.to.plot <- unique(protein.traces.with.complex$complex_id)

pdf('corum_complexes.pdf', width=7.5, height=6)
idx <- 1
for (complex.id in complex.ids.to.plot) {
    cat(paste(sprintf('plotting complex %s', complex.id),
              sprintf('(%d / %d)', idx, length(complex.ids.to.plot)),
              sep='\t'),
        '\n')
    traces <- protein.traces.with.complex[complex_id == complex.id, ]
    complex.name <- traces[1, complex_name]
    p <- ggplot(traces) +
         geom_line(aes(x=sec, y=intensity, color=protein_id)) +
         # theme(legend.position='none') +
         ggtitle(sprintf('%s | %s', complex.id, complex.name)) +
         scale_x_continuous(breaks=round(seq(0, max(traces$sec), by=10), 1))
    annotation.label <- paste(
        traces[1, n_proteins_in_complex], '/',
        traces[1, n_proteins_in_complete_complex],
        sep=''
    )
    max.intensity <- traces[, max(intensity)]
    max.sec <- traces[, max(sec)]
    p <- p + annotate('text', x=0.95*max.sec, y=0.96*max.intensity,
                      label=annotation.label, size=12)
    print(p)
    idx <- idx + 1
}
dev.off()

# Write the ids and names of the complexes to a TSV file
complex.overview.table <-
    unique(protein.traces.with.complex[,
           list(complex_id, complex_name, n_proteins_in_complex,
                n_proteins_in_complete_complex)])
complex.overview.table[['apexes (fully_observed)']] <- ''
complex.overview.table[['apexes (partially_observed)']] <- ''
write.table(complex.overview.table, 'corum_complexes_overview.tsv',
            sep='\t', row.names=F)
