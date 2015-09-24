######################################################################
#                              PREAMBLE                              #
######################################################################
# Set the current directory
setwd('/Users/robin/Desktop/MA_SideTask/')

# Load all required packages, if a packge is not installed already, install it
required.packages <- c('ggplot2', 'data.table')
for (p in required.packages) {
    import.ok <- require(p, char=TRUE)
    if (!import.ok) {
        chooseCRANmirror(ind=81)  # choose ETH by index
        install.package(p)
        require(p, char=TRUE)
    }
}

######################################################################
#                            DATA IMPORT                             #
######################################################################
# Read the pairwise correlation data (NCI60)
nci60.raw <- fread('data/NCI60/correlations.tsv')
setnames(nci60.raw, c('prot1', 'prot2', 'r1', 'r2'))
# Remove second column with correlation coefficients since it is identical
# to the first one.
nci60 <- subset(nci60.raw, select=-r2)
setnames(nci60, 'r1', 'r')

# Get a association list protein id <-> corum complex id
complex.assoc <- fread('data/corum_complex_protein_assoc.tsv')
# Drop the unecessary name column
complex.assoc <- subset(complex.assoc, select=-complex_name)
complex.assoc[, unique(protein_id), by=complex_id]

######################################################################
#           PRODUE INTER-COMPLEX PROTEIN TRACE CORRELATION           #
######################################################################
# Read in a long list file of peptide traces, the dataframe has the columns:
# ['protein_id', 'peptide_id', 'sec', 'peptide_intensity']
peptide.traces <- fread('data/HEK293_peptide_traces_long.tsv')

# Produce the protein trace by summing up the individual traces (group by both
# protein id and sec number, then sum up the intensity within the groups).
protein.traces <- peptide.traces[, sum(peptide_intensity),
                                 by=list(protein_id, sec)]
setnames(protein.traces, 'V1', 'intensity')


# Assign to each protein trace observation the complex id to which
# the observation's protein belongs.
setkey(complex.assoc, protein_id)  # set key for later merge operation
setkey(protein.traces, protein_id)
protein.tr.annot <- merge(protein.traces, complex.assoc, allow.cartesian=T)

all.obs.ids <- unique(protein.tr.annot$complex_id)  # the observed proteins
setkey(complex.assoc, complex_id)
setkey(protein.tr.annot, complex_id)

# Check for each complex if all of its proteins are present
is.full.complex <- sapply(all.obs.ids, function(cid) {
    all.prots <- complex.assoc[cid, protein_id]
    obs.prots <- protein.tr.annot[cid, protein_id]
    is.full.compl <- setequal(all.prots, obs.prots)
    is.full.compl
})
is.complete.flag <- data.table(complex_id=all.obs.ids, is_full=is.full.complex,
                               key='complex_id')
# Merge the above computed flags into the trace data frame. Then we can
# drop all the rows that have is_full == FALSE (i.e., the row/intensity observation
# belongs to a protein that belongs to a complex which isn't fully observed).
protein.tr.annot.flagged <- protein.tr.annot[is.complete.flag]
prot.tr <- protein.tr.annot.flagged[is_full == TRUE]
prot.tr[, is_full := NULL]  # delete the flag inplace since we don't need it anymore

prot.tr.wide <- dcast.data.table(prot.tr, protein_id + complex_id ~ sec,
                                 value.var='intensity')

setkey(prot.tr.wide, complex_id)
complex.ids <- unique(prot.tr.wide$complex_id)
avg.corr <- sapply(complex.ids, function(cid) {
    # Get all protein traces that belong to the complex with id `cid`
    prot.traces <- prot.tr.wide[cid]
    # Drop non-numeric columns and transpose the dataframe, producing
    # a dataframe where rows are sec numbers and columns are protein ids
    trace.mat <- t(subset(prot.traces, select=-c(protein_id, complex_id)))
    # Compute pairwise correlations between traces and take the mean
    mean(cor(trace.mat))
})
complex.avg.corr <- data.table(complex_id=complex.ids, avg_corr=avg.corr)

######################################################################
#                     PRODUCE NCI60 CORRELATIONS                     #
######################################################################

