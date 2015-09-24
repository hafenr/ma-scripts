import pandas as pd
import os
from corum_data import complex_assoc, is_complete_complex

path_to_script = os.path.dirname(os.path.abspath(__file__))

# Read the pairwise correlation data (NCI60), takes some time...
nci60_raw = pd.read_csv(
    os.path.join(path_to_script, 'data/NCI60/correlations.tsv'), sep='\t',
    header=None, names=['protein_id_1', 'protein_id_2', 'r1', 'r2']
)

# Drop second correlation column (axis=1) and rename it since pairwise correlation
# is the same regardless of argument order.
nci60 = nci60_raw.drop('r2', axis=1).rename(columns={'r1': 'r'})

# Since the imported pairwise correlations are computed for each protein pair
# (even if they aren't in the same complex), we need to remove all those
# pairwise correlations that are between proteins that aren't in the same complex.
complex_id_assoc = complex_assoc[['protein_id', 'complex_id']]
# First annotate the first protein with its complex membership
prot_annot_ = pd.merge(nci60, complex_id_assoc,
                       left_on='protein_id_1', right_on='protein_id').\
                 rename(columns={'complex_id': 'complex_id_1'})

# Then also annotate the second protein with its complex membership
prot_annot = pd.merge(prot_annot_, complex_id_assoc,
                      left_on='protein_id_2', right_on='protein_id').\
                rename(columns={'complex_id': 'complex_id_2'})
# Only select those rows that are pairwise correlations between proteins
# that belong to the same complex.
belongs_to_same_complex = prot_annot.complex_id_1 == prot_annot.complex_id_2
prot_annot_same_complex = prot_annot.loc[belongs_to_same_complex]
# Drop some of the columns introduced by the merge and clean up the column names
prot_annot_same_complex_ = \
    prot_annot_same_complex[['protein_id_1', 'protein_id_2', 'complex_id_1', 'r']].\
    rename(columns={'complex_id_1': 'complex_id'})

def filter_complete_complexes(df):
    protein_ids = set(df.protein_id_1).union(set(df.protein_id_2))
    return is_complete_complex(df.complex_id.iloc[0], protein_ids)

prot_annot_same_complex = \
    prot_annot_same_complex_.groupby('complex_id').\
                             filter(filter_complete_complexes)

# Finally, group all rows together that have the same complex id and compute
# the mean.
complex_corr_nci60 = prot_annot_same_complex.groupby('complex_id').mean().r
