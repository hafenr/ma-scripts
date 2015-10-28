require(data.table)


# peptr.file <- '~/Dev/MAScripts/data/peptide_traces_with_complex_id.tsv'
prot.conc.estimates.file <- '~/Dev/MAScripts/data/HEK293_protein_conc_estimates.csv'

peptraces <- fread('~/Desktop/cprophet/input_peptides_cprophet.csv', header=T)
peptide.sums <- rowSums(subset(peptraces, select=-c(protein_id, peptide_id)))
peptraces$peptide_sum <- peptide.sums
peptraces <- subset(peptraces, select=c(protein_id, peptide_sum))
prot.intensity <- peptraces[, sum(peptide_sum), by=protein_id]
setnames(prot.intensity, 'V1', 'total_intensity')


# Sum peptide traces together to produce the protein traces
# peptr <- fread(peptr.file)
# prot.intensity <- setnames(peptr[, sum(intensity), by=list(protein_id)],
#                            'V1', 'total_intensity')
prot.with.conc <- fread(prot.conc.estimates.file)

# Merge the DT holding the estimated protein concentrations and the DT holding
# the summed up intensities across the SEC dimension.
# The new DT shouldn't contain protein ids that aren't present in either of the DTs
# that are to be merged (nomatch=0).
setkey(prot.intensity, protein_id)
setkey(prot.with.conc, protein_id)
merged <- prot.intensity[prot.with.conc, nomatch=0]

# Fit a line to the logged values
fit <- lm(log(protein_concentration) ~ log(total_intensity), merged)
prots.imputed <- subset(prot.intensity,
                        !(prot.intensity$protein_id %in% prot.with.conc$protein_id))
# Inpute the concentration values for these proteins (log concentrations are converted
# back to their original unit).
prots.imputed[, protein_concentration := exp(predict(fit, prots.imputed))]

# Plot how this stuff looks like
with(merged, plot(log(total_intensity), log(protein_concentration)))

p <- smooth.spline(log(merged$total_intensity),
                   log(merged$protein_concentration))
fitted.vals <- predict(p, log(merged$total_intensity))
points(fitted.vals$x, fitted.vals$y, col='blue')
with(prots.imputed, points(log(total_intensity), log(protein_concentration),
                           col='green', pch=4))

abline(fit, col='red', lwd=3)

# For each of the proteins with an estimated conc. there is also an entry indicating
# the proteins molecular weight. So for the proteins with imputed concentrations, we also
# need to get those molecular weight values.
# List from: http://www.uniprot.org/uploadlists/ by pasting output from
# paste0(prots.imputed$protein_id, collapse=' ') into the text field
proteins.without.mws <-  paste0(prots.imputed$protein_id, collapse=' ')
mw.file <- '~/Dev/MAScripts/data/HEK293_uniprot_missing_molweights.tsv'
molweights <- fread(mw.file)[, list(Entry, Mass)]
molweights <- setnames(molweights, c('protein_id', 'protein_mw'))
# Since uniprot saves their numbers with commas in them they are read as strings and
# we need to convert them back.
molweights[, protein_mw := as.numeric(gsub(',', '', protein_mw))]
# Add the molecular weight to the other proteins with imputed values
setkey(molweights, protein_id)
setkey(prots.imputed, protein_id)
prots.imputed <- prots.imputed[molweights]
# Delete unused column
prots.imputed <- prots.imputed[, total_intensity := NULL]

# Combine the proteins that had a concentration estimates with those that now have an
# imputed value and write to file.
all.proteins <- rbind(prot.with.conc, prots.imputed)
write.table(all.proteins, '~/Dev/MAScripts/data/HEK293_all_proteins_with_concentration.csv',
            row.names=F, sep=',')

