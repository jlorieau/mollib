"""
MolLib settings file.

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-04T20:29:39-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-05T09:57:58-05:00
   @License:            Copyright 2016
"""

# Ramachandran parameters

# The optimal Ramachandran angles for a helix. The threshold is the range (+/-)
# That an angle has to be within the optimal values to be considered a helix.
dihedral_helix_angles = (-60, -45)
dihedral_helix_angles_threshold = (50, 50)


# Hydrogenation parameters

# The optimal bond length (in Angstroms) for various bonds.
bond_length = {'N-H': 1.023,    # L. Yao, et al. JACS 130, 16518-20 (2008)
               'C-H': 1.117,    # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'CA-HA': 1.117,  # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'C-H3': 1.106,   # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'C-C': 1.517,    # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               }


# Hydrogen bond parameters

# The cutoff distance (in Angstroms) between donor and acceptor atoms to be
# considered a valid hydrogen bond
hbond_cutoff = 2.5

# The default threshold angle range (+/-) for hydrogen bonds
hbond_angle_threshold = 30

# alpha-helix hydrogen bond target CO-HN angle and threshold
hbond_a_helix_angle = 149.
hbond_a_helix_angle_threshold = hbond_angle_threshold

# 310-helix hydrogen bond target CO-HN angle and threshold
hbond_310_helix_angle = 114.
hbond_310_helix_angle_threshold = hbond_angle_threshold

# pi-helix hydrogen bond target CO-HN angle and threshold
hbond_pi_helix_angle = 149.
hbond_pi_helix_angle_threshold = hbond_angle_threshold

# beta-sheet hydrogen bond target CO-HN angle and threshold
hbond_beta_sheet_angle = 155.
hbond_beta_sheet_angle_threshold = hbond_angle_threshold
