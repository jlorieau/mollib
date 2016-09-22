#: Element selector string for the hydrogen bond donor atom1. The OR
#: operator, '|', is supported.
#:
#: Hydrogen bond dipoles are defined as:
#:
#:     donor2--donor1 .... acceptor1--acceptor2
donor1_elements = 'H|D'

#: Element selector string for the hydrogen bond donor atom2. The OR
#: operator, '|', is supported.
donor2_elements = 'N|15N|O'

#: Element selector string for the hydrogen bond acceptor atom1. The OR
#: operator, '|', is supported.
acceptor1_elements = 'O'

#: Element selector string for the hydrogen bond acceptor atom2. The OR
#: operator, '|', is supported.
acceptor2_elements = 'C|13C'


#: The cutoff distance ranges (in A) between atoms to be considered a hydrogen
#: bond.
#:
#: 'd1a1': The distance between the donor1 and acceptor1 atoms
hbond_distance_cutoff = {'d1a1': (1.8, 2.8),
                         }

#: The cutoff angle ranges (in degrees) between atoms to be considered a
#: hydrogen bond.
hbond_angle_cutoff = {'theta': (110., 180.),
                       'phi': (-180., 180.)
                      }

# Phi torsion angle range (in degrees) for helices. (Generously allowed)
helix_phi = (-170., 0.)

# Phi torsion angle range (in degrees) for helices. (Generously allowed)
helix_psi = (-100., 55.)

# Phi torsion angle range (in degrees) for beta-sheet. (Generously allowed)
beta_phi = (-200.,-25.)

# Psi torsion angle range (in degrees) for beta-sheet. (Generously allowed)
beta_psi = (45., 220.)


#: Phi torsion angle ranges (in degrees) for residue i+1 in turns
beta_turn_i1_phi = {'turnI':   (-110., -10.),
                    'turnII':  (-110., -10.),
                    'turnIp':  (10., 110.),
                    'turnIIp': (10., 110.),
                    }

#: Psi torsion angle ranges (in degrees) for residue i+1 in turns
beta_turn_i1_psi = {'turnI':   (-80., 20.),
                    'turnII':  (70., 170.),
                    'turnIp':  (-20., 80.),
                    'turnIIp': (-170., -70.)
                    }

#: Phi torsion angle ranges (in degrees) for residue i+2 in turns
beta_turn_i2_phi =  {'turnI': (-140., -40.),
                     'turnII': (30., 140.),
                     'turnIp': (40., 140.),
                     'turnIIp': (-130., -30.)
                     }

#: Psi torsion angle ranges (in degrees) for residue i+2 in turns
beta_turn_i2_psi = {'turnI':   (-50., 50.),
                    'turnII':  (-50., 50.),
                    'turnIp':  (-50., 50.),
                    'turnIIp': (-50., 50.)
                    }

#: Major classification name for backbone-backbone (bb-bb) amide hydrogen
#: bonds
major_bb_bb_amide = 'bb-bb amide'

#: Major classification name for backbone-sidechain (sc-bb) amide hydrogen
#: bonds
major_bb_sc_amide = 'bb-sc amide'

#: Major classification name for sidechain-backbone (sc-bb) amide hydrogen
#: bonds
major_sc_bb_amide = 'sc-bb amide'

#: Major classification name for sidechain-sidechain (sc-sc) amide hydrogen
#: bonds
major_sc_sc_amide = 'sc-sc amide'

#: Major classification name for backbone-backbone (bb-bb) aliphatic hydrogen
#: bonds
major_bb_bb_aliphatic = 'bb-bb aliph.'

#: Major classification name for backbone-sidechain (bb-sc) aliphatic hydrogen
#: bonds
major_bb_sc_aliphatic = 'bb-sc aliph.'

#: Major classification name for sidechain-backbone (sc-bb) aliphatic hydrogen
#: bonds
major_sc_bb_aliphatic = 'sc-bb aliph.'

#: Major classification name for backbone-backbone (sc-sc) aliphatic hydrogen
#: bonds
major_sc_sc_aliphatic = 'sc-sc aliph.'

#: Major classification name for backbone-backbone (bb-sc) hydroxyl hydrogen
#: bonds
major_bb_sc_hydroxyl = 'bb-sc hydroxyl'

#: Major classification name for sidechain-backbone (sc-bb) hydroxyl hydrogen
#: bonds
major_sc_bb_hydroxyl = 'sc-bb hydroxyl'

#: Major classification name for sidechain-sidechain (sc-sc) hydroxyl hydrogen
#: bonds
major_sc_sc_hydroxyl = 'sc-sc hydroxyl'

#: Minor classification name for type I turns
minor_beta_turnI = "type I turn"

#: Minor classification name for type II turns
minor_beta_turnII = "type II turn"

#: Minor classification name for type I' turns
minor_beta_turnIp = "type I' turn"

#: Minor classification name for type II' turns
minor_beta_turnIIp = "type II' turn"

#: Minor classification name for sheets
minor_beta = 'sheet'

#: Minor classification name for anti-parallel beta-sheets
minor_beta_anti = 'sheet, anti-parallel'

#: Minor classification name for parallel beta-sheets
minor_beta_par = 'sheet, parallel'

#: Minor classification name for 310-helices
minor_310 = '310-helix'

#: Minor classification name for alpha-helices
minor_alpha = 'alpha-helix'

#: Minor classification name for pi-helices
minor_pi = 'pi-helix'

#: Minor classification name for isolated hydrogen bonds
minor_isolated = 'isolated'
