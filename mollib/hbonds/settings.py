donor1_elements = 'H|D'
donor2_elements = 'N|15N|O'
"""Element strings for the donor elements. The '|' represents the 'OR'
operator.

  Hydrogen bond dipoles are defined as:

    donor2--donor1 .... acceptor1--acceptor2
"""

acceptor1_elements = 'O'
acceptor2_elements = 'C|13C'
"""Element strings for the acceptor elements. The '|' represents the 'OR'
operator.
"""

hbond_distance_cutoff = {'d1a1': (1.8, 2.5),
                         }
"""The cutoff distance ranges (in A) between atoms to be considered a hydrogen
bond.

'd1a1':
    The distance between the donor1 and acceptor1 atoms
"""

hbond_angle_cutoff = {'theta': (110., 180.),
                       'phi': (-180., 180.)
                      }
"""The cutoff angle ranges (in deg) between atoms to be considered a
hydrogen bond.

'd2a1a2':
    The angle between donor2--acceptor1--acceptor2
"""