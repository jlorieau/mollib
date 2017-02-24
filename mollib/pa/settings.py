#: The gyromagnetic ratios (in rad T^-1 s^-1) of common nuclei
gamma = {'H': 267.513E6,
         'C': 67.262E6,
         'N': -27.116E6}

default_predicted_rdcs = (('N', 'H'),
                          ('CA', 'HA'),
                          ('CA', 'HA2'),
                          ('CA', 'HA3'),
                          ('N', 'C'),
                          ('CA', 'C'),
                          )

default_predicted_racs = {'H': {'dxx': -5.0,  # ppm
                                }
                          }
