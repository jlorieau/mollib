"""
MolLib default settings
"""
# Author: Justin L Lorieau
# Copyright: 2016

# TODO: move plugin settings to their own plugin submodule

# Ramachandran parameters

dihedral_helix_angles = (-60, -45)
"The default optimal Ramachandran angles for a helix (in degrees)."

dihedral_helix_angles_threshold = (30, 30)
"""The default Ramanchandran angle threshold (+/- threshold) to be considered a
helix (in degrees)"""

# Hydrogenation parameters

bond_length = {'N-H': 1.023,    # L. Yao, et al. JACS 130, 16518-20 (2008)
               'C-H': 1.117,    # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'CA-HA': 1.117,  # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'C-H3': 1.106,   # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               'C-C': 1.517,    # M. Ottiger, et al. JACS 121, 4690-4695 (1999)
               }
"The default optimal length of standard bonds (in Angstroms) for biomolecules"

# Hydrogen bond parameters

hbond_cutoff = 2.5
"The default acceptor-donor atom distance to be considered a hydrogen bond."

hbond_amide_cutoff = hbond_cutoff
"""The default acceptor-donor atom distance to be considered an *amide*
hydrogen bond."""

hbond_aliphatic_cutoff = 2.8
"""The default acceptor-donor atom distance to be considered and *aliphatic*
hydrogen bond."""

hbond_angle_threshold = 30.
"""The default angle threshold (in degrees)to be considered a valid hydrogen
bond (+/- threshold)."""

hbond_a_helix_coh_angle = 149.
"The default hydrogen bond C-O-HN angle for an *alpha-helix*."

hbond_a_helix_angle_threshold = hbond_angle_threshold
"""The default hydrogen bond angle threshold from the ideal angle of an
*alpha-helix*."""

hbond_310_helix_coh_angle = 114.
"The default hydrogen bond C-O-HN angle for a *310-helix*."

hbond_310_helix_angle_threshold = hbond_angle_threshold
"""The default hydrogen bond angle threshold from the ideal angle of a
*310-helix*."""

hbond_pi_helix_coh_angle = 149.
"The default hydrogen bond C-O-HN angle for a *pi-helix*."

hbond_pi_helix_angle_threshold = hbond_angle_threshold
"""The default hydrogen bond angle threshold from the ideal angle of a
*pi-helix*."""

hbond_beta_sheet_coh_angle = 155.
"The default hydrogen bond C-O-HN angle for a *beta-sheet*."

hbond_beta_sheet_angle_threshold = hbond_angle_threshold
"""The default hydrogen bond angle threshold from the ideal angle of a
*beta-sheet*."""