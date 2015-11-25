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

complex.features <- fread('~/Desktop/cc_out/sec_proteins_with_dscore_filtered.csv')
complex.features.fil <- complex.features[d_score > 2, ]


protein.traces.wo.complex <- subset(protein.traces.long, select=-c(complex_name, complex_id))
setkey(protein.traces.wo.complex)
protein.traces.wo.complex <- unique(protein.traces.wo.complex)

complex.features.fil <- 
    complex.features.fil[order(complex.features.fil$d_score, decreasing=T), ]
for (i in 1:nrow(complex.features.fil)) {
    feat <- complex.features.fil[i, ]
    cat('DECOY: ', feat$decoy == 1, '\n')
    cat('D_SCORE: ', feat$d_score, '\n')
    subunits <- strsplit(feat$aggr_Fragment_Annotation, ';')[[1]]
    trace.data <- protein.traces.wo.complex[protein_id %in% subunits, ]
    feat$max_intensity <- max(trace.data$intensity)
    p <- ggplot(aes(x=sec, y=intensity, color=protein_id),
           data=trace.data) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept=feat$leftWidth) +
        geom_vline(xintercept=feat$rightWidth)
    print(p)
    readline(prompt = "Pause. Press <Enter> to continue...")
}
