#: The gyromagnetic ratios (in rad T^-1 s^-1) of common nuclei
gamma = {'H': 267.513E6,
         'C': 67.262E6,
         'N': -27.116E6}

#: Calculate dipolar couplings from bond lengths and gyromagnetic ratios. If
#: this is False, the values in the default_predicted_rdcs are used. However,
#: if a value is not specified in the default_predicted_rdcs, the value will
#: still be calculated from bond lengths.
calculate_from_bonds = False

#: The default dipolar couplings calculated.
#: The static dipolar coupling constants listed are for common bonds in
#: proteins. Values are calculated based on average calculated DCCs for 2MJB
#: and scaled to an HN DCC of -11472Hz. This is equal to the
#: dipolar reduced anisotropy.
#: These values are used if ``calculate_from_bonds`` is False.
default_predicted_rdcs = {'N-H': 10823., #11472.,  # 1.02 A
                          'NE1-HE1': 10823.,
                          'CA-HA': -22300.,
                          'C-N+1': 1115.,
                          'C-CA': -1880.,
                          'CA-CB': -1880.,
                          }


#: The default RACS tensor values of various nuclei in proteins.
#: order - the xx/yy/zz component order of the tensor.
#:
#:          1. The first component is colinear with the atom -- ref_atom1
#:             vector
#:
#:          2. The second component is orthogonal to the
#:             atom -- ref_atom1 -- ref_atom2 plane
#:
#:          3. The third component is orthogonal to the second two components.
#:
#:          4. Tensor rotations conducted with the Z-Y-X convention for the
#:             alpha, beta and gamma angles, respectively.
#:
#:  See the following reference for tensor conventions and values: Cornilescu
#:  et al. JACS 2000, 122, 10143.

default_predicted_racs = {
    # The Carbonyl CSA. This convention places the x-axis (s11) near the,
    # C-N bond and the z-axis (s33) orthogonal to the O-C-N plane.
    'C': {'delta': -86.53 * 1.03,  # ppm (Reduced anisotropy)
          'eta': 0.63,
          'alpha': 40.,            # degrees
          'beta': 0.,              # degrees
          'gamma': 0.,             # degrees
          'ref_atom1': 'N+1',
          'ref_atom2': 'O',
          'order': 'xzy',
          },

    # The Nitrogen CSA. This convention places the z-axis (s11) near the H-N
    # bond, and the y-axis (s22) orthogonal to the H-N-CA plane.
    'N': {'delta': 108.53 * 0.993,  # ppm (Reduced anisotropy)
          'eta': 0.16,
          'alpha': 0.,              # degrees
          'beta': -20.,             # degrees
          'gamma': 0.,              # degrees
          'ref_atom1': 'H',
          'ref_atom2': 'CA',
          'order': 'zyx',
          },

    # The amide H CSA. This convention places the z-axis (s33) near the H-N
    # bond, and the x-axis (s11) orthogonal to the H-N-C plane.
    'H': {'delta': -5.93,           # ppm (Reduced anisotropy)
          'eta': 1.00,
          'alpha': 0.,              # degrees
          'beta': 0.,               # degrees
          'gamma': -7.,              # degrees
          'ref_atom1': 'N',
          'ref_atom2': 'C-1',
          'order': 'zxy',
          },
    }

#: Default errors in absolute values
default_error = {'N-H':   0.2,  # Hz
                 'CA-HA': 0.5,  # Hz
                 'N-C':   0.5,  # Hz
                 'CA-C':  0.5,  # HZ
                 'NE-HE': 0.2,  # Hz
                 'C-C':   0.5,  # Hz (bond type)
                 'C-H':   1.0,  # Hz (bond type)
                 'C':     1.,  # ppb
                 'N':     1.,  # ppb
                 'H':     5.,  # ppb
                 }

#: urls to download mr files
mr_urls = ('https://files.rcsb.org/download/',)

#: Xplor-NIH incorporates C-H methyl RDCs by projecting them on the associated
#: carbon-carbon bond. If the data uses this convention, then set this setting
#: to True. This convention, however, occludes the use of C-C RDCs, like the
#: CA-CB coupling of an ALA, and it is not set by default.
project_methyls = False

#: The order parameter to use for methyl groups
methyl_order_parameter = 1.0

#: Settings to enable individual fixers by default
enable_signfixer = True
enable_outlierfixer = False
enable_nhscalefixer = False

#: The alpha-critical values for the Grubbs test of outliers. These are used to
#: identify warning or bad data points in terms of their deviations with
#: respect to the best-fit SVD data.
alpha_warning = 0.05
alpha_bad = 0.01


