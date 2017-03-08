#: The gyromagnetic ratios (in rad T^-1 s^-1) of common nuclei
gamma = {'H': 267.513E6,
         'C': 67.262E6,
         'N': -27.116E6}

#: Calculate dipolar couplings from bond lengths and gyromagnetic ratios. If
#: this is False, the values in the default_predicted_rdcs are used.
calculate_from_bonds = True

#: The default dipolar couplings and static dipolar coupling constants for
#: common bonds in proteins. Values are calculated based on average calculated
#: DCCs for 2MJB and scaled to an HN DCC of -11472Hz. These values are used
#: if ``calculate_from_bonds`` is False.
default_predicted_rdcs = {('N', 'H'): 11472.,  # 1.02 A
                          ('CA', 'HA'): -21500.,
                          ('CA', 'HA2'): -21500.,
                          ('CA', 'HA3'): -21500.,
                          ('N', 'C-1'): 1115.,
                          ('CA', 'C'): -1880.
                          }

default_predicted_racs = {'C': {'delta': -86.53,  # ppm
                                'eta': 0.73,
                                'alpha': -12.5,     # degrees
                                'beta': 0.,
                                'ref_atom1': 'O',
                                'ref_atom2': 'N+1',
                                }
                          }

optimize_to_interactions = (('N', 'H'),
                            ('CA', 'HA'),
                            )

#: The alpha-critical values to identify warning or bad data points in terms
#: of their deviations with respect to the best-fit SVD data.
alpha_warning = 0.05
alpha_bad = 0.01
