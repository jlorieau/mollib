#: The gyromagnetic ratios (in rad T^-1 s^-1) of common nuclei
gamma = {'H': 267.513E6,
         'C': 67.262E6,
         'N': -27.116E6}

#: Calculate dipolar couplings from bond lengths and gyromagnetic ratios. If
#: this is False, the values in the default_predicted_rdcs are used.
calculate_from_bonds = False

#: The default dipolar couplings calculated.
#  The static dipolar coupling constants listed are for common bonds in
#: proteins. Values are calculated based on average calculated DCCs for 2MJB
#: and scaled to an HN DCC of -11472Hz. This is equal to the
#: dipolar reduced anisotropy.
#: These values are used if ``calculate_from_bonds`` is False.
default_predicted_rdcs = {'N-H': 10823., #11472.,  # 1.02 A
                          'NE1-HE1': 10823.,
                          'CA-HA': -22300.,
                          'N-C-1': 1115.,
                          'CA-C': -1880.
                          }


#: order - the xx/yy/zz component order of the tensor.
#:
#:          1. The first component is colinear with the atom -- ref_atom1
#:             vector
#:          2. The second component is orthogonal to the
#:             atom -- ref_atom1 -- ref_atom2 plane
#:          3. The third component is orthogonal to the second two components.
default_predicted_racs = {
    'C': {'delta': -86.53 * 1.03,  # ppm (Reduced anisotropy)
          'eta': 0.63,
          'alpha': 40.,     # degrees
          'beta': 0.,
          'gamma': 0.,
          'ref_atom1': 'N+1',
          'ref_atom2': 'O',
          'order': 'xzy',
          },
    'N': {'delta': 108.53 * 0.993,  # ppm (Reduced anisotropy)
          'eta': 0.16,
          'alpha': 0.,     # degrees
          'beta': -20.,
          'gamma': 0.,
          'ref_atom1': 'H',
          'ref_atom2': 'CA',
          'order': 'zyx',
          },
    }

#: Default errors in absolute values
default_error = {'N-H':   0.2,  # Hz
                 'CA-HA': 0.5,  # Hz
                 'N-C':   0.5,  # Hz
                 'CA-C':  0.5,  # Hz
                 'NE-HE': 0.2,  # Hz
                 'C':     1.,  # ppb
                 'N':     1.,  # ppb
                 }

optimize_to_interactions = (('N', 'H'),
                            ('CA', 'HA'),
                            )

#: The alpha-critical values to identify warning or bad data points in terms
#: of their deviations with respect to the best-fit SVD data.
alpha_warning = 0.05
alpha_bad = 0.01

#: Settings to enable individual fixers by default
enable_signfixer = True
enable_outlierfixer = False
