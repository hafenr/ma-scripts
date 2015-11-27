require(devtools)
require(ggplot2)
require(hash)
# MACode git commit: f80708bd55d2d7724573a4d56599d078ce27e736
devtools::load_all('~/Dev/MACode')


annotations.1.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/151109_SEC-SWATH_CORUM_manual_annotation_mhe.tsv'
annotations.2.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/COMPLEXES_4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED_ANNOTATED.tsv'

manual.annotations.1.raw <- readManualAnnotationFile(annotations.1.raw)
manual.annotations.2.raw <- readManualAnnotationFile(annotations.2.raw)

complex.completeness.info <- manual.annotations.1.raw$complexes
complete.complex.ids <- complex.completeness.info[
    n_proteins_in_complex == n_proteins_in_complete_complex,
    complex_id
]

manual.annotations.1 <-
    createManualComplexAnnotations(manual.annotations.1.raw$annotations, 
                                   'apexes_partially_observed') 
manual.annotations.2 <-
    createManualComplexAnnotations(manual.annotations.2.raw$annotations,
                                   'apexes_partially_observed') 

manual.annotations <- mergeManualComplexAnnotations(manual.annotations.1,
                                                    manual.annotations.2)

manual.annotations.complete <-
    manual.annotations.full[complex_id %in% complete.complex.ids, ]
