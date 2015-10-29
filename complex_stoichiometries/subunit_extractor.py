"""
subunit_extractor
-----------------

A short script to extract all biological units (assemblies) from a directory
of PDB files that were downloaded from RCSB PDB.

Use the script like this:

    $ python subunit_extractor path/to/download/dir \
            -o . \
            --corum path/to/corum/file.tsv

"""

import os.path as p
import pandas as pd
from collections import defaultdict
import string
from subprocess import call
import os
import os.path as p
import re
import argparse
import logging


logger = logging.getLogger(__name__)

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT, level=logging.DEBUG)


class Chain:
    """Container class for PDB chains."""
    def __init__(self, letter, trans_nr, protein, seq_start, seq_end):
        self.letter = letter
        self.trans_nr = trans_nr
        self.protein = protein
        self.seq_start = seq_start
        self.seq_end = seq_end

    def __str__(self):
        return '%s:%s[%d:%d]' % (
            self.letter, self.protein, self.seq_start, self.seq_end)


class BiolUnit:
    """
    Container class for biological units as they are described in PDB files
    under section REMARK 350

    """
    def __init__(self, pdb_id, nr, chains):
        """Summary line.

        Parameters
        ----------
        pdb_id : string
            The pdb id to which this biological unit belongs.
        nr : number
            The number of the biological unit as reported in REMARK 350.
        chains : List[Chain]
            A list of chains. A chain in this sense is a chain as reported
            by the PDB file, but which has a specific spatial orientation.

        Returns
        -------
        BiolUnit

        """
        self.pdb_id = pdb_id
        self.nr = nr
        self.chains = chains

        all_chains_identified = \
            all([ch.protein is not None for ch in self.chains])
        if all_chains_identified:
            proteins = set([ch.protein for ch in self.chains])
            n_prots = len(proteins)

            # Stoichiometries should look like A2B2 etc.
            prot_to_letter = \
                dict(zip(proteins, string.ascii_uppercase[:n_prots]))

            # TODO: Some complexes will show up as having a homomer stoichiometry A2, even though
            # they are classified as heteromers. I assume this happends because all the DBREFs
            # have the same protein id in them. They differ however in their seq_start/seq_end attribute.
            # So I assume the heteromer is built from different subsequences of a protein seq.
            # One would need to treat them specially but I think its ok to just treat them as homomers.
            # One could also compute their molecular weight based on the subsequences instead of the full protein sequences
            # (this might be a good idea in general, since many chains correspond to a subsequence of a uniprot sequence).
            self.stoichiometry = defaultdict(int)
            for ch in self.chains:
                self.stoichiometry[prot_to_letter[ch.protein]] += 1

            # Stoichiometry describing string as found on the PDB website
            stoich_string = \
                ''.join([''.join(map(str, s))
                        for s in self.stoichiometry.items()])
            self.stoich_string = stoich_string
        else:
            self.proteins = None
            self.stoichiometry = None
            self.stoich_string = ''

    def __str__(self):
        return '%s - %d \t chains: %s \t stoich: %s' % (
            self.pdb_id, self.nr, ', '.join(map(str, self.chains)),
            self.stoich_string)

    @property
    def is_homomer(self):
        return len(self.stoichiometry.keys()) == 1

    def as_dataframe(self):
        """Convert this biological unit to a pandas dataframe where each row
        corresponds to a chain that belongs to this unit."""
        chains = []
        for ch in self.chains:
            chains.append({
                'pdb_id': self.pdb_id,
                'biol_unit_nr': self.nr,
                'letter': ch.letter,
                'trans_nr': ch.trans_nr,
                'protein_id': ch.protein,
                'seq_start': ch.seq_start,
                'seq_end': ch.seq_end,
                'complex_stoichiometry': self.stoich_string
            })
        return pd.DataFrame(chains)


def extract_biological_units(pdb_id, text, dbrefs):
    """Read the section on biological assemblies form a PDB file and output
    a list of objects of class BiolUnit.
    Note that biological units that contain chains that have no DBREF entry
    which points to Uniprot are ignored.

    Parameters
    ----------
    pdb_id : string
        The pdb id that corresponds to the file being analyzed.
    text : string
        The contents of the file as a single string.
    dbrefs : Dict[string, Dict]
        A dictionary that maps chain letters to another dictionary that
        holds the information contained in the DBREF section.

    Returns
    -------
    List[BiolUnit]

    """
    # Be careful with newline characters. Each linebreak in the multiline string is treated
    # as a \n in the resulting regexp. Regexes for targeting repeats of lines
    # therefore need to have an optional \n appended to them.
    regex = re.compile(r"""
REMARK 350 BIOMOLECULE: (?P<biomolecule_nr>\d)\s*
(REMARK 350 \w.*\n?)*?
REMARK 350 APPLY THE FOLLOWING TO CHAINS: (?P<apply_to_chains>.*)
(?P<trans_lines>REMARK 350   BIOMT1.*
(.|\n)*?
REMARK 350   BIOMT3.*)
(REMARK \d\d\d)\s*$""",
        re.MULTILINE)  # ^$ match start/end w.r.t to \n instead of whole string.

    def match_to_biol_unit(m):
        nr = int(m.group('biomolecule_nr'))
        trans_lines = m.group('trans_lines').split('\n')
        apply_to_chains = map(str.strip, m.group('apply_to_chains').split(','))
        n_trans = len(trans_lines) / 3
        # Not all chains have a row in the DBREF section, i.e. there is no
        # Uniprot ID for this chain => return None
        if not all([ch in dbrefs for ch in apply_to_chains]):
            return None
        else:
            chains = []
            # Create a Chain object for each chain with a specific spatial orientation.
            # This means that chains are duplicated according to the transformations that
            # should be applied to them to get the biological unit.
            for chain_letter in apply_to_chains:
                for i in range(1, n_trans + 1):  # duplicate
                    protein = dbrefs[chain_letter]['protein']
                    seq_start = int(dbrefs[chain_letter]['seq_start'])
                    seq_end = int(dbrefs[chain_letter]['seq_end'])
                    ch = Chain(chain_letter, i, protein, seq_start, seq_end)
                    chains.append(ch)
            bunit = BiolUnit(pdb_id, nr, chains)
            print bunit
            return bunit

    matches = list(regex.finditer(text))
    biolunits = []
    for m in matches:
        bunit = match_to_biol_unit(m)
        # If bunit is None, not all chains have a DBREF entry.
        if bunit is not None:
            biolunits.append(bunit)
    return biolunits


def get_dbrefs(text):
    """
    Extract the DBREF portion of the PDB file and return it as a dictionary
    that maps each chain letter to another dictionary that contains the
    respective fields.
    Only DBREFs that map to uniprot are considered.

    """

    regex = re.compile(
        r"^DBREF\s+(?P<pdb_id>\w+)\s+(?P<chain_id>\w+)\s+(?P<seq_begin>\d+)\s+(?P<insert_begin>\d+)\s+(?P<db_id>\w+)\s+(?P<db_acc>\w+)\s+(?P<seq_db_code>\w+)\s+(?P<seq_start>\d+)\s+(?P<seq_end>\d+).*$",
        re.MULTILINE
    )
    info = {}
    matches = regex.finditer(text)
    for m in matches:
        chain_id = m.group('chain_id')
        if m.group('db_id') == 'UNP':
            info[chain_id] = {
                'protein': m.group('db_acc'),
                'seq_start': m.group('seq_start'),
                'seq_end': m.group('seq_end')
            }
    return info


def get_all_chains_from_files_in_dir(input_dir):
    """Load all pdb files (potentially zipped) from directory a directory
    and extract the biological units.

    Parameters
    ----------
    input_dir : string
        The path to the directory holding the downloaded pdb files.

    Returns
    -------
    pd.DataFrame
        A dataframe where each row corresponds to a chain of some biological
        unit.

    """
    data_dir = input_dir
    pdb_files = [p.join(data_dir, f) for f in os.listdir(data_dir)]
    complete_biounits = []

    for fname in pdb_files:
        logger.info('current pdb file %s' % fname)
        if fname.endswith('.gz') and not p.exists(p.splitext(fname)[0]):
            call(['gzip', '-d', fname])
            logger.info('zipped pdb file extracted' % fname)
            fname = p.splitext(fname)[0]
        with open(fname, 'r') as f:
            text = f.read()
            pdb_id = p.basename(p.splitext(fname)[0])
            dbrefs = get_dbrefs(text)
            bunits = extract_biological_units(pdb_id, text, dbrefs)
            logger.info('found %d complete units' % len(bunits))
            complete_biounits += bunits
    all_chains_df = pd.concat(
        [bunit.as_dataframe() for bunit in complete_biounits])

    n_pdb_ids = len(all_chains_df.pdb_id.drop_duplicates())
    logger.info(
        'found a total of %d complete biological units belonging '
        'to %d pdb entries'
        % (len(complete_biounits), n_pdb_ids))

    return all_chains_df


def merge_with_corum_complexes(chains_df, corum_complexes):
    """Merge a df of chains extracted from PDB files with a list
    of corum complexes.

    Parameters
    ----------
    chains_df : pd.DataFrame
        A dataframe where each row is a chain belonging to some PDB
        biological unit. The rows must have the columns:
        - pdb_id                 # the pdb entry id
        - pdb_biol_unit          # a number indicating the biological unit
        - protein_id             # a uniprot identifier
    corum_complexes : pd.DataFrame
        A dataframe of CORUM complex <-> Uniprot ID associations.
        The DF should have at least the following columns:
        - complex_id             # CORUM identifier
        - complex_name           # some string describing the complex
        - protein_id             # uniprot identifier
        - complex_stoichiometry  # a string such as A2B2

    Returns
    -------
    pd.DataFrame

    """
    # Add a new column containing the unique protein ids as a frozen set
    # to both the corum complexes and the pdb biological units.
    # Since frozensets are hashable, we can later use this column
    # to merge both dataframes.
    def add_unique_proteins(df):
        df['unique_proteins'] = [frozenset(df.protein_id)] * len(df)
        return df
    corum_complexes_with_unique = \
        corum_complexes.groupby('complex_id').apply(add_unique_proteins)
    # Drop duplicates so that rows are now complexes instead of complex-protein
    # associations.
    corum_complexes_with_unique = \
        corum_complexes_with_unique[['complex_id', 'complex_name', 'unique_proteins']].\
        drop_duplicates()

    def add_stoichs(df):
        df['unique_proteins'] = [frozenset(df.protein_id)] * len(df)
        # Sort the proteins first so that the two protein lists A,B and B,A
        # both will be added as 'A,B', 'A,'B'
        proten_list = sorted(df.protein_id)
        df['protein_stoichiometry'] = [','.join(proten_list)] * len(df)
        return df
    pdb_chains_with_unique = chains_df.groupby(['pdb_id', 'biol_unit_nr']).apply(add_stoichs)
    # Identical units of a complex are thrown out at this point!
    pdb_complexes_with_unique = \
        pdb_chains_with_unique[['pdb_id', 'protein_stoichiometry', 'complex_stoichiometry', 'unique_proteins']].\
        drop_duplicates()


    complexes_with_stoich = pd.merge(corum_complexes_with_unique, pdb_complexes_with_unique)
    logger.info(
        'was able to find one or more stoichiometries for %d corum complexes'
        % len(complexes_with_stoich.complex_id.drop_duplicates()))
    logger.info(
        'average number of stoichs per complex: %f'
        % complexes_with_stoich.groupby('complex_id').\
                                complex_id.count().mean())

    return complexes_with_stoich.drop('unique_proteins', axis=1)


def main(input_dir, output_dir, corum_complexes_file):
    logger.info('start program')
    chains_df = get_all_chains_from_files_in_dir(input_dir)

    chains_output_file = 'pdb_chains.tsv'
    chains_df.to_csv(chains_output_file, index=False, sep='\t')
    logger.info('saved all chains to %s' % chains_output_file)

    if corum_complexes_file:
        logger.info('merging with corum')
        corum_complexes = pd.read_csv(corum_complexes_file, sep='\t')
        merged_complexes = merge_with_corum_complexes(
            chains_df, corum_complexes)
        corum_output_file = 'corum_complexes_with_stoichs.tsv'

        merged_complexes.to_csv(corum_output_file, sep='\t', index=False)
        logger.info('saved all corum complexes to %s' % corum_output_file)
    else:
        logger.info('no corum file supplied, won\'t merge')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract biological units form pdb files.')
    parser.add_argument(
        'input_dir', help=(
            'input directory holding pdb files as downloaded by the '
            'RCSB PDB java downloader (possibly zipped)'))
    parser.add_argument(
        '-o', '--out',
        help=('directory where to store the TSV file holding the individual '
              'chains  and the file holding the CORUM complexes'))
    parser.add_argument(
        '-c', '--corum',
        help='a file with columns complex_id, protein_id, and complex_name')
    args = parser.parse_args()

    main(args.input_dir, args.out, args.corum)
