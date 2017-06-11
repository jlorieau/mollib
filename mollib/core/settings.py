#: The default pH of new molecules
default_pH = 7.0

#: The default pKs of ionizable amino-acids. [Ref]_
#:
#: Some amino-acids have degenerate ionizeable atoms; these are listed and
#: separated by '-' characters. The different pKs for each ionization is
#: listed in the subsequent items in the tuple.
#:
#:  .. [Ref] G. R. Grimsley, J. M. Scholtz, C. N. Pace, Protein Sci. 18, 247-51
#:           (2009).
pKs = {'ASP': {'OD1-OD2': (-1.0, 3.5)},
       'GLU': {'OE1-OE2': (-1.0, 4.2)},
       'HIS': {'ND1-NE2': (6.6, 14.0)},
       'CYS': {'SG': (6.8,)},
       'TYR': {'OH': (10.3,)},
       'LYS': {'NZ': (10.5, 14.0, 14.0)},
       'last': {'O-OXT': (-1.0, 3.3)},
       'first': {'N': (7.7, 14.0, 14.0),}
     }

#: Path for the datasets
dataset_path = 'data/'

#: Path for model input molecule identifiers
model_molecule_identifiers = ('high_res.txt',
                              #'high_res_short.txt',
                              )

#: Path for the Ramachandran statistics datasets
ramachandran_dataset_path = 'data/ramachandranstatistics/'

#: Path for the Hbond statistics datasets
hbond_dataset_path = 'data/hbondstatistics/'

#: The cutoff Energy (in kT) to report good, warning and bad energies.
energy_cutoff_good = 3.4  # Within 96.6% of observed values
energy_cutoff_warning = 5.4  # Within 99.5% of observed values
energy_cutoff_bad = 20.0  # Only observed < 0.5% of the time.

#: urls to download PDB files
pdb_urls = ('https://files.rcsb.org/download/',)

#: Only load the first model, by default
pdb_first_model = True
