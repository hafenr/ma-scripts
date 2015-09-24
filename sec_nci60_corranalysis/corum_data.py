import pandas as pd
import os

path_to_script = os.path.dirname(os.path.abspath(__file__))

# Get a association list protein id <-> corum complex id
complex_assoc = pd.read_csv(
    os.path.join(path_to_script, 'data/corum_complex_protein_assoc.tsv'),
    sep='\t'
).drop('complex_name', axis=1)
# Given the above list, construct a hash map that links every complex id to the
# set of protein identifies that consitute this complex.
complex_ids = complex_assoc.complex_id.drop_duplicates()
protein_sets = complex_assoc.groupby('complex_id').\
                             apply(lambda df: set(df.protein_id)).tolist()
complexes = dict(zip(complex_ids, protein_sets))

def is_complete_complex(complex_id, protein_ids):
    return complexes[complex_id] == set(protein_ids)
