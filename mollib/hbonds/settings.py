donor1_elements = 'H|D'
donor2_elements = 'N|15N|O'

acceptor1_elements = 'O'
acceptor2_elements = 'C|13C'

hbond_distance_cutoff = {'d2a2': (2.0, 2.5),
                        }
"""The cutoff distance ranges (in A) between atoms to be considered a hydrogen
bond.

'd2a2':
    The distance between the donor2 and acceptor2 atoms
"""
hbond_angle_cutoff_ = {'d2a1a2': (100., 160.)
                      }
"""The cutoff angle ranges (in deg) between atoms to be considered a
hydrogen bond.

'd2a1a2':
    The angle between donor2--acceptor1--acceptor2
"""