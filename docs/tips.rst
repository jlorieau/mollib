Tips and Tricks
###############

Shell Aliases and Functions
---------------------------

The following aliases are helpful in minimizing typing at the command line.
These can be added to the user shell script (ex: ``.bashrc``).

.. code:: shell

    # mollib: backbone dihedral angles
    alias mlr='mollib measure --rama -i'

    # mollib: hydrogen bonds
    alias mlh='mollib hbonds --hydrogenate -i'

    # mollib: partial alignment from the PDB
    function mlp() { mollib pa -i $1 -a $1 ${@:2};}

    # mollib: download a PDB file
    function mlget() { mollib process -i $1 --save; gzip -d $1.pdb.gz;}
