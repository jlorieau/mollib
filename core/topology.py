
# Bonding Topology of heavy atoms. This is a dict of a dict of sets
topology = {'PRO': {'N':   {'C-1', 'CA', 'CD'},
                    'CA':  {'N', 'C', 'CB', 'HA'},
                    'CB':  {'CA', 'CG', 'HB2', 'HB3'},
                    'CG':  {'CB', 'CD', 'HG2', 'HG3'},
                    'CD':  {'CG', 'N', 'HD2', 'HD3'},
                    'C':   {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'GLY': {'N':   {'C-1', 'CA', 'H'},
                    'CA':  {'N', 'C', 'HA2', 'HA3'},
                    'C':   {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ALA': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB1', 'HB2', 'HB3'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ARG': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'CD': {'CG', 'HD2', 'HD3', 'NE'},
                    'NE': {'CD', 'HE', 'CZ'},
                    'CZ': {'NE', 'NH1', 'NH2',},
                    'NH1': {'CZ', 'HH11', 'HH12'},
                    'NH2': {'CZ', 'HH21', 'HH22'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ASN': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'OD1', 'ND2'},
                    'OD1': {'CG'},
                    'ND2': {'CG', 'HD21', 'HD22'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ASP': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'OD1', 'OD2'},
                    'OD1': {'CG','HD1'},
                    'OD2': {'CG','HD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'CYS': {'N': {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'SG', 'HB2', 'HB3'},
                    'SG': {'CB', 'HG'},
                    'C': {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'GLN': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'CD': {'CG', 'OE1', 'NE2'},
                    'OE1': {'CD'},
                    'NE2': {'CG', 'HE21', 'HE22'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'GLU': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'CD': {'CG', 'OE1', 'OE2'},
                    'OE1': {'CD', 'HE1'},
                    'OE2': {'CD', 'HE2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'HIS': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'ND1', 'CD2'},
                    'ND1': {'CG', 'CE1', 'HD1'},
                    'CE1': {'ND1', 'NE2', 'HE1'},
                    'NE2': {'CE1', 'CD2', 'HE2'},
                    'CD2': {'CG', 'NE2', 'HD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ILE': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB', 'CG1', 'CG2'},
                    'CG1': {'CB', 'HG12', 'HG13', 'CD1'},
                    'CD1': {'CG1', 'HD1', 'HD2', 'HD3'},
                    'CG2': {'CB', 'HG1', 'HG2', 'HG3'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'LEU': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'CD1', 'CD2', 'HG'},
                    'CD1': {'CG', 'HD11', 'HD12', 'HD13'},
                    'CD2': {'CG', 'HD21', 'HD22', 'HD23'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'LYS': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'CD': {'CG', 'HD2', 'HD3', 'CE'},
                    'CE': {'CD', 'HE2', 'HE3', 'NZ'},
                    'NZ': {'CE', 'HZ1', 'HZ2', 'HZ3'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'MET': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'HG2', 'HG3', 'SD'},
                    'SD': {'CG', 'CE'},
                    'CE': {'SD', 'HE1', 'HE2', 'HE3'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'PHE': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'CD1', 'CD2'},
                    'CD1': {'CG', 'HD1', 'CE1'},
                    'CE1': {'CD1', 'HE1', 'CZ'},
                    'CZ': {'CE1', 'HZ', 'CE2'},
                    'CE2': {'CZ', 'HE2', 'CD2'},
                    'CD2': {'CE2', 'CD2', 'CG'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'SER': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'OG'},
                    'OG': {'CB', 'HG'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'THR': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB', 'CG2', 'OG1'},
                    'OG1': {'CB', 'HG1'},
                    'CG2': {'CB', 'HG21', 'HG22', 'HG23'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'TRP': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA',  'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'CD1', 'CD2'},
                    'CD1': {'CG', 'HD1', 'NE1'},
                    'NE1': {'CD1', 'HE1', 'CE2'},
                    'CE2': {'NE1', 'CD2', 'CZ2'},
                    'CZ2': {'CE2', 'HZ2', 'CH2'},
                    'CH2': {'CZ2', 'HH2', 'CZ3'},
                    'CZ3': {'CH2', 'HZ3', 'CE3'},
                    'CE3': {'CZ3', 'HE3', 'CD2'},
                    'CD2': {'CE3', 'CE2', 'CG'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'TYR': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'CG': {'CB', 'CD1', 'CD2'},
                    'CD1': {'CG', 'HD1', 'CE1'},
                    'CE1': {'CD1', 'HE1', 'CZ'},
                    'CZ': {'CE1', 'CE2', 'OH'},
                    'OH': {'CZ', 'HH'},
                    'CE2': {'CZ', 'HE2', 'CD2'},
                    'CD2': {'CE2', 'HD2', 'CG'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'VAL': {'N':  {'C-1', 'CA', 'H'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'CB': {'CA', 'HB', 'CG1', 'CG2'},
                    'CG1': {'CB', 'HG11', 'HG12', 'HG13'},
                    'CG2': {'CB', 'HG21', 'HG22', 'HG23'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            }