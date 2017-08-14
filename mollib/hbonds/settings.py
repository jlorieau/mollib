# Hydrogen Bond Definitions
# -------------------------

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
hbond_distance_cutoff = {'d1a1': (1.5, 2.5),
                         }

#: The cutoff angle ranges (in degrees) between atoms to be considered a
#: hydrogen bond.
hbond_angle_cutoff = {'theta': (105., 180.),
                       'phi': (-180., 180.)
                      }


# Dihedral Angle Ranges
# ---------------------

# Phi torsion angle range (in degrees) for helices. (Generously allowed)
helix_phi = (-170., 0.)

# Phi torsion angle range (in degrees) for helices. (Generously allowed)
helix_psi = (-100., 55.)

# Phi torsion angle range (in degrees) for beta-sheet. (Generously allowed)
beta_phi = (-200., -25.)

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


# Classification names
# --------------------

#: Classification type name for backbone-backbone (bb-bb) amide hydrogen
#: bonds
type_bb_bb_amide = 'bb-bb amide'

#: Classification type name for backbone-sidechain (sc-bb) amide hydrogen
#: bonds
type_bb_sc_amide = 'bb-sc amide'

#: Classification type name for sidechain-backbone (sc-bb) amide hydrogen
#: bonds
type_sc_bb_amide = 'sc-bb amide'

#: Classification type name for sidechain-sidechain (sc-sc) amide hydrogen
#: bonds
type_sc_sc_amide = 'sc-sc amide'

#: Classification type name for backbone-backbone (bb-bb) aliphatic hydrogen
#: bonds
type_bb_bb_aliphatic = 'bb-bb aliph.'

#: Classification type name for backbone-sidechain (bb-sc) aliphatic hydrogen
#: bonds
type_bb_sc_aliphatic = 'bb-sc aliph.'

#: Classification type name for sidechain-backbone (sc-bb) aliphatic hydrogen
#: bonds
type_sc_bb_aliphatic = 'sc-bb aliph.'

#: Classification type name for backbone-backbone (sc-sc) aliphatic hydrogen
#: bonds
type_sc_sc_aliphatic = 'sc-sc aliph.'

#: Classification type name for backbone-backbone (bb-sc) hydroxyl hydrogen
#: bonds
type_bb_sc_hydroxyl = 'bb-sc hydroxyl'

#: Classification type name for sidechain-backbone (sc-bb) hydroxyl hydrogen
#: bonds
type_sc_bb_hydroxyl = 'sc-bb hydroxyl'

#: Classification type name for sidechain-sidechain (sc-sc) hydroxyl hydrogen
#: bonds
type_sc_sc_hydroxyl = 'sc-sc hydroxyl'

#: Major classification name for type I turns
major_beta_turnI = "type I turn"

#: Major classification name for type II turns
major_beta_turnII = "type II turn"

#: Major classification name for type I' turns
major_beta_turnIp = "type I' turn"

#: Major classification name for type II' turns
major_beta_turnIIp = "type II' turn"

#: Major classification name for sheets
major_beta = 'sheet'

#: Major classification name for anti-parallel beta-sheets
major_beta_anti = 'sheet, anti-parallel'

#: Major classification name for parallel beta-sheets
major_beta_par = 'sheet, parallel'

#: Major classification name for 310-helices
major_310 = '310-helix'

#: Major classification name for alpha-helices
major_alpha = 'alpha-helix'

#: Major classification name for pi-helices
major_pi = 'pi-helix'

#: Major classification name for isolated hydrogen bonds
major_isolated = 'isolated'

#: Minor classification name for the N-terminal residues
minor_N = 'N-term'

#: Minor classification name for the C-terminal residues of alpha-helices
minor_C = 'C-term'

#: Minor classification name for glycines
minor_gly = 'Gly'


# Assign Blocks
# -------------

#: Assign blocks of residues, including gaps
assign_blocks = True

#: Assign alpha-helix blocks
assign_blocks_alpha = True

#: Assign beta-strand blocks, including gaps
assign_blocks_beta = True

#: Assign 310-helix blocks, including gaps
assign_blocks_310 = True

#: assign_blocks: Overwrite classification assignments, if an assignment has
#: already been made
assign_blocks_overwrite = False


#: assign_blocks: Check the previous and subsequent residues in a contiguous
#: alpha-helix assignment
assign_blocks_alpha_extend_termini = False

#: assign_blocks: Label the 'minor' classification of the specified number of
# N- or C-terminal residues as 'N-term' or 'C-term'
assign_blocks_alpha_label_N_term = 1
assign_blocks_alpha_label_C_term = 1

#: assign_blocks: Only fill gaps in alpha-helices that are as large as the
#: following number
assign_blocks_alpha_gap_tolerance = 1


#: assign_blocks: Check the previous and subsequent residues in a contiguous
#: beta-sheet assignment
assign_blocks_beta_extend_termini = True

#: assign_blocks: Label the 'minor' classification of the specified number of
# N- or C-terminal residues as 'N-term' or 'C-term'
assign_blocks_beta_label_N_term = 1
assign_blocks_beta_label_C_term = 1

#: assign_blocks: Only fill gaps in beta-sheets that are as large as the
#: following number
assign_blocks_beta_gap_tolerance = 2


#: assign_blocks: Check the previous and subsequent residues in a contiguous
#: 310 assignment
assign_blocks_310_extend_termini = False

#: assign_blocks: Label the 'minor' classification of the specified number of
#: N- or C-terminal residues as 'N-term' or 'C-term'
assign_blocks_310_label_N_term = 0
assign_blocks_310_label_C_term = 0

#: assign_blocks: Only fill gaps in 310-helices that are as large as the
#: following number
assign_blocks_310_gap_tolerance = 3


# Tables
# ------

#: Render HBond tables with detailed information
hbond_table_detailed = False

#: Sort HBond tables by the hydrogen bond type
hbond_table_sort_type = False

#: Render the Rama table with detailed information, including minor assignments
rama_table_detailed = True
