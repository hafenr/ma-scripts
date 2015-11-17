require(devtools)
# MACode git commit: f80708bd55d2d7724573a4d56599d078ce27e736
devtools::load_all('~/Dev/MACode')


annotations.1.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/151109_SEC-SWATH_CORUM_manual_annotation_mhe.tsv'
annotations.2.raw <-
    '~/Dev/MAScripts/manual_complex_annotation/COMPLEXES_4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED_ANNOTATED.tsv'

manual.annotations <- mergeManualComplexAnnotations(annotations.1.raw, annotations.2.raw,
                                   'apexes_fully_observed')

# sec_complexes output form cprophet
detected.features <-
    fread('~/Dev/cprophet/run-674/iteration-0000/sec_complexes.tsv')

# Only look at the manual features were every measured 
# protein was also part of the peak (doesn't have to be
# a complete protein group).
true.features <- manual.annotations[apex_type == 'apexes_fully_observed']

assessed.feats <- assessComplexFeatures(true.features, detected.features, feature.vicinity.tol=5)
