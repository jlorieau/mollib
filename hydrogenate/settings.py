
bond_length = {'N-H': 1.023,    # L. Yao, et al. JACS 130, 16518-20 (2008)
               'C-H': 1.117,    # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'CA-HA': 1.117,  # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'C-H3': 1.106,   # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'C-C': 1.517,    # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'N-H2': 1.01,    # PK Sawinski, et al. Crys Growth & Design
                                # 13, 1730 (2013)
               }
"""The default optimal length of standard bonds (in Angstroms) for biomolecules
"""

amide_atom_name = 'HN'  # FIXME: Test changing this name
"""The default name to use for creating new amide protons."""

default_pH = 7.0
"""The default pH to use when adding hydrogens to molecules.
"""

pKs = {'ASP': 3.5,
       'GLU': 4.2,
       'HIS': 6.6,
       'CYS': 6.8,
       'TYR': 10.3,
       'LYS': 10.5,
       'C-term': 3.3,
       'N-term': 7.7}
"""The default pKs of ionizable amino-acids. [Ref]_


  .. [Ref]: G. R. Grimsley, J. M. Scholtz, C. N. Pace,
            Protein Sci. 18, 247-51 (2009).
"""
