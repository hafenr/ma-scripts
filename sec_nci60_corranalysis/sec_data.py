import pandas as pd
import numpy as np
import os
from corum_data import complex_assoc, is_complete_complex

path_to_script = os.path.dirname(os.path.abspath(__file__))
os.chdir(path_to_script)

# Read in a long list file of peptide traces, the dataframe has the columns:
# [u'protein_id', u'peptide_id', u'sec', u'peptide_intensity']
peptide_traces = pd.read_csv('data/HEK293_peptide_traces_long.tsv', sep='\t').\
                    rename(columns={'peptide_intensity': 'intensity'})

# Sum up all peptide traces to produce the protein trace
protein_traces_ = peptide_traces.groupby(['protein_id', 'sec']).\
                                sum().\
                                reset_index()

# Get a list of all protein ids that are were identified in the sample
protein_ids_in_sample = protein_traces_[['protein_id']].drop_duplicates()
# Annotate the protein ids with the complex id they could be involved
protein_ids_in_sample_annot = \
    protein_ids_in_sample.merge(complex_assoc[['complex_id', 'protein_id']])
# Get a list of all protein ids that belong to complete complexes.
# (Note that this only means that the protein was measured AT SOME POINT, not
# necessarily that the protein eluted when the complex eluted from the SEC
# column).
def check_if_complex_is_complete(df):
    complex_id = df.complex_id.iloc[0]
    return is_complete_complex(complex_id, df.protein_id)
protein_ids_of_complete_complexes = \
    protein_ids_in_sample_annot.groupby('complex_id').\
                                filter(check_if_complex_is_complete).\
                                protein_id.drop_duplicates()

# Convert the long list data frame to a wide list data frame, similar to a
# matrix.
# Rows now correspond to vectors of intensities. The index (row labels) are the
# protein ids.
protein_traces = protein_traces_.loc[
    protein_traces_.protein_id.isin(protein_ids_of_complete_complexes)
]

protein_traces_mat = protein_traces.pivot(
    index='protein_id', values='intensity', columns='sec'
)

# Annotate each row (protein trace) with a complex label that is taken from
# the dataframe complex assox
protein_traces_mat_annot_ = protein_traces_mat.merge(
    complex_assoc, right_on='protein_id', left_index=True
)

# Due to the merge above, the complex id and protein id were put back into the
# dataframe as columns. Let's remove them and use them as hierarchical row
# labels. E.g.:
#                    sec1  sec2  sec3 .... sec80
# complex_1 prot_1a   0     1     2          0
#           prot_1b   0     5     4          0
# complex_2 prot_2a   ...
#           prot_2b   ...
protmat_with_monos = \
    protein_traces_mat_annot_.set_index(['complex_id', 'protein_id'])

# Remove all complexes that only have one protein
protmat = protmat_with_monos.groupby(level='complex_id').\
                             filter(lambda df: len(df) != 1)

def compute_correlation(intensity_df):
    # Get the values as a 2d numpy array
    mat = intensity_df.values
    # Compute the correlations between all row vectors (= protein traces of
    # a specific complex)
    corrmat = np.corrcoef(mat)
    n_proteins = len(intensity_df)
    # Extract the upper triangular section of the matrix
    # (with an offset of one above the diagonal)
    triu_idx = np.triu_indices(n_proteins, 1)
    # Calculate the average correlation
    avg_r = corrmat[triu_idx].mean()
    return avg_r
# Produce a series of average correlation values, the index (value labels) are
# the CORUM complex ids.
complex_corr_sec = protmat.groupby(level='complex_id').apply(compute_correlation)
