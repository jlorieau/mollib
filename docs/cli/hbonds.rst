.. _hbonds_command:

Hbonds Command
==============
The ``hbonds`` command detects and reports hydrogen bonds in molecules.
The ``hbonds`` command has the following features:

- It identifies hydrogen bonds between backbones, sidechains, subunits and
  ligands.
- It has a flexible interface to identify a wide-range of hydrogen bond types
  including *amide*, *aliphatic*, *hydroxyl* and the specification of
  arbritrary electric dipole types.
- It classifies backbone-backbone hydrogen bonds based on backbone torsion
  angles.

    - Helical stretches require contiguous residues with helical backbone
      torsion angles. This prevents the misclassification of 310-helices and
      beta turns as well as isolated i+3, i+4 and i+5 hydrogen bonds.

Usage
-----

    .. include:: output/cli_hbonds_help.html


    .. note:: If the molecule does not have hydrogens, this command will need
              to be run with the ``--hydrogenate`` parameter. See the
              :ref:`hydrogenate <hydrogenate>` option.

Arguments
---------

    ``--aliphatic``
        (Optional) Include aliphatic hydrogen bonds in the results. The
        acceptor-donoratom cutoff distances are elongated to 3.0 A, and carbon
        atoms are allowed in hydrogen bond donor dipoles.

    ``--detailed``
        (Optional) Present a detailed report on the geometries of hydrogen
        bonds. Does not include classification information.

    ``--sort-type``
        (Optional) Sort the hydrogen bonds by classification type, like
        beta-sheet or alpha-helix.


Examples
--------

    .. include:: output/cli_hbonds_i_2KXA.html
