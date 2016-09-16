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

helix_phi = (-170., 0.)
helix_psi = (-90., 45.)
"""Generously allowed torsion angle ranges (in degrees) for helices.
"""

beta_phi = (-180.,-25.)
beta_psi = (-50., 180.)
"""Generously allowed torsion angle ranges (in degrees) for beta sheets."""

beta_turn_i1_phi = {'turnI':   (-110., -10.),
                    'turnII':  (-110., -10.),
                    'turnIp':  (10., 110.),
                    'turnIIp': (10., 110.),
                    }

beta_turn_i1_psi = {'turnI':   (-80., 20.),
                    'turnII':  (70., 170.),
                    'turnIp':  (-20., 80.),
                    'turnIIp': (-170., -70.)
                    }

beta_turn_i2_phi =  {'turnI': (-140., -40.),
                     'turnII': (30., 140.),
                     'turnIp': (40., 140.),
                     'turnIIp': (-130., -30.)
                     }

beta_turn_i2_psi = {'turnI':   (-50., 50.),
                    'turnII':  (-50., 50.),
                    'turnIp':  (-50., 50.),
                    'turnIIp': (-50., 50.)
                    }
"""Torsion angle ranges (in degrees) for beta turns.
"""
