
# Bonding Topology of heavy atoms. This is a dict of a dict of sets
topology = {'PRO': {'N':   {'C-1', 'CA', 'CD'},
                    'CA':  {'N', 'C', 'CB', 'HA'},
                    'HA':  {'CA'},
                    'CB':  {'CA', 'CG', 'HB2', 'HB3'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG':  {'CB', 'CD', 'HG2', 'HG3'},
                    'HG2': {'CG'},
                    'HG3': {'CG'},
                    'CD':  {'CG', 'N', 'HD2', 'HD3'},
                    'HD2': {'CD'},
                    'HD3': {'CD'},
                    'C':   {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'GLY': {'N':   {'C-1', 'CA', 'H'},
                    'H':   {'N'},
                    'CA':  {'N', 'C', 'HA2', 'HA3'},
                    'HA2': {'CA'},
                    'HA3': {'CA'},
                    'C':   {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ALA': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB1', 'HB2', 'HB3'},
                    'HB1': {'CB'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ARG': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'HG2': {'CG'},
                    'HG3': {'CG'},
                    'CD': {'CG', 'HD2', 'HD3', 'NE'},
                    'HD2': {'CD'},
                    'HD3': {'CD'},
                    'NE': {'CD', 'HE', 'CZ'},
                    'HE': {'NE'},
                    'CZ': {'NE', 'NH1', 'NH2',},
                    'NH1': {'CZ', 'HH11', 'HH12'},
                    'HH11': {'NH1'},
                    'HH12': {'NH1'},
                    'NH2':  {'CZ', 'HH21', 'HH22'},
                    'HH21': {'NH2'},
                    'HH22': {'NH2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ASN': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG':  {'CB', 'OD1', 'ND2'},
                    'OD1': {'CG'},
                    'ND2': {'CG', 'HD21', 'HD22'},
                    'HD21': {'ND2'},
                    'HD22': {'ND2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ASP': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'OD1', 'OD2'},
                    'OD1': {'CG','HD1'},
                    'HD1': {'OD1'},
                    'OD2': {'CG','HD2'},
                    'HD2': {'OD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'CYS': {'N': {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CB'},
                    'CB': {'CA', 'SG', 'HB2', 'HB3'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'SG': {'CB', 'HG'},
                    'HG': {'SG'},
                    'C': {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'GLN': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'HG2': {'CG'},
                    'HG3': {'CG'},
                    'CD': {'CG', 'OE1', 'NE2'},
                    'OE1': {'CD'},
                    'NE2': {'CD', 'HE21', 'HE22'},
                    'HE21': {'NE2'},
                    'HE22': {'NE2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'GLU': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'HG2': {'CG'},
                    'HG3': {'CG'},
                    'CD': {'CG', 'OE1', 'OE2'},
                    'OE1': {'CD', 'HE1'},
                    'HE1': {'OE1'},
                    'OE2': {'CD', 'HE2'},
                    'HE2': {'OE2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'HIS': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'ND1', 'CD2'},
                    'ND1': {'CG', 'CE1', 'HD1'},
                    'HD1': {'ND1'},
                    'CE1': {'ND1', 'NE2', 'HE1'},
                    'HE1': {'CE1'},
                    'NE2': {'CE1', 'CD2', 'HE2'},
                    'HE2': {'NE2'},
                    'CD2': {'CG', 'NE2', 'HD2'},
                    'HD2': {'CD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'ILE': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB', 'CG1', 'CG2'},
                    'HB': {'CB'},
                    'CG1': {'CB', 'HG12', 'HG13', 'CD1'},
                    'HG12': {'CG1'},
                    'HG13': {'CG1'},
                    'CD1': {'CG1', 'HD11', 'HD12', 'HD13'},
                    'HD11': {'CD1'},
                    'HD12': {'CD1'},
                    'HD13': {'CD1'},
                    'CG2': {'CB', 'HG21', 'HG22', 'HG23'},
                    'HG21': {'CG2'},
                    'HG22': {'CG2'},
                    'HG23': {'CG2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'LEU': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'CD1', 'CD2', 'HG'},
                    'HG': {'CG'},
                    'CD1': {'CG', 'HD11', 'HD12', 'HD13'},
                    'HD11': {'CD1'},
                    'HD12': {'CD1'},
                    'HD13': {'CD1'},
                    'CD2': {'CG', 'HD21', 'HD22', 'HD23'},
                    'HD21': {'CD2'},
                    'HD22': {'CD2'},
                    'HD23': {'CD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'LYS': {'N':  {'C-1', 'CA', 'H'},
                    'H':  {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'HG2', 'HG3', 'CD'},
                    'HG2': {'CG'},
                    'HG3': {'CG'},
                    'CD': {'CG', 'HD2', 'HD3', 'CE'},
                    'HD2': {'CD'},
                    'HD3': {'CD'},
                    'CE': {'CD', 'HE2', 'HE3', 'NZ'},
                    'HE2': {'CE'},
                    'HE3': {'CE'},
                    'NZ': {'CE', 'HZ1', 'HZ2', 'HZ3'},
                    'HZ1': {'NZ'},
                    'HZ2': {'NZ'},
                    'HZ3': {'NZ'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'MET': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'HG2', 'HG3', 'SD'},
                    'HG2': {'CG'},
                    'HG3': {'CG'},
                    'SD': {'CG', 'CE'},
                    'CE': {'SD', 'HE1', 'HE2', 'HE3'},
                    'HE1': {'CE'},
                    'HE2': {'CE'},
                    'HE3': {'CE'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'PHE': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'CD1', 'CD2'},
                    'CD1': {'CG', 'HD1', 'CE1'},
                    'HD1': {'CD1'},
                    'CE1': {'CD1', 'HE1', 'CZ'},
                    'HE1': {'CE1'},
                    'CZ': {'CE1', 'HZ', 'CE2'},
                    'HZ': {'CZ'},
                    'CE2': {'CZ', 'HE2', 'CD2'},
                    'HE2': {'CE2'},
                    'CD2': {'CE2', 'HD2', 'CG'},
                    'HD2': {'CD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'SER': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'OG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'OG': {'CB', 'HG'},
                    'HG': {'OG'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'THR': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB', 'CG2', 'OG1'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'OG1': {'CB', 'HG1'},
                    'HG1': {'HG1'},
                    'CG2': {'CB', 'HG21', 'HG22', 'HG23'},
                    'HG21': {'CG2'},
                    'HG22': {'CG2'},
                    'HG23': {'CG2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'TRP': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA',  'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'CD1', 'CD2'},
                    'CD1': {'CG', 'HD1', 'NE1'},
                    'HD1': {'CD1'},
                    'NE1': {'CD1', 'HE1', 'CE2'},
                    'HE1': {'NE1'},
                    'CE2': {'NE1', 'CD2', 'CZ2'},
                    'CZ2': {'CE2', 'HZ2', 'CH2'},
                    'HZ2': {'CZ2'},
                    'CH2': {'CZ2', 'HH2', 'CZ3'},
                    'HH2': {'CH2'},
                    'CZ3': {'CH2', 'HZ3', 'CE3'},
                    'HZ3': {'CZ3'},
                    'CE3': {'CZ3', 'HE3', 'CD2'},
                    'HE3': {'CE3'},
                    'CD2': {'CE3', 'CE2', 'CG'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'TYR': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB2', 'HB3', 'CG'},
                    'HB2': {'CB'},
                    'HB3': {'CB'},
                    'CG': {'CB', 'CD1', 'CD2'},
                    'CD1': {'CG', 'HD1', 'CE1'},
                    'HD1': {'CD1'},
                    'CE1': {'CD1', 'HE1', 'CZ'},
                    'HE1': {'CE1'},
                    'CZ': {'CE1', 'CE2', 'OH'},
                    'OH': {'CZ', 'HH'},
                    'HH': {'OH'},
                    'CE2': {'CZ', 'HE2', 'CD2'},
                    'HE2': {'CE2'},
                    'CD2': {'CE2', 'HD2', 'CG'},
                    'HD2': {'CD2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            'VAL': {'N':  {'C-1', 'CA', 'H'},
                    'H': {'N'},
                    'CA': {'N', 'C', 'CB', 'HA'},
                    'HA': {'CA'},
                    'CB': {'CA', 'HB', 'CG1', 'CG2'},
                    'HB': {'CB'},
                    'CG1': {'CB', 'HG11', 'HG12', 'HG13'},
                    'HG11': {'CG1'},
                    'HG12': {'CG1'},
                    'HG12': {'CG1'},
                    'CG2': {'CB', 'HG21', 'HG22', 'HG23'},
                    'HG21': {'CG2'},
                    'HG22': {'CG2'},
                    'HG22': {'CG2'},
                    'C':  {'CA', 'O', 'N+1'},
                    'O':   {'C'},
                    },
            }