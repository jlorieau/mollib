# Bonding Topology. This is a dict of a dict of sets
topology = {'PRO': {'N':   {'C-1', 'CA', 'CD'},
                            'CA':  {'N', 'C', 'CB', 'HA'},
                            'CB':  {'CA', 'CG', 'HB2', 'HB3'},
                            'CG':  {'CB', 'CD', 'HG2', 'HG3'},
                            'CD':  {'CG', 'N', 'HD2', 'HD3'},
                            'C':   {'CA', 'N+1'}
                            },
                    'GLY': {'N':   {'C-1', 'CA', 'HN'},
                            'CA':  {'N', 'C', 'HA'},
                            'HA2': {'CA'},
                            'HA3': {'CA'},
                            'C':   {'CA', 'N+1'},
                            },
                    'ALA': {'N':  {'C-1', 'CA', 'HN'},
                            'CA': {'N', 'C', 'CB', 'HA'},
                            'CB': {'CA', 'HB1', 'HB2', 'HB3'},
                            'C':  {'CA', 'N+1'},
                            },
                    'ASP': {'OD1': {'CG','HD1'},
                            'OD2': {'CG','HD2'},
                            },

                    }