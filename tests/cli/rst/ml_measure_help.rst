.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span> ml measure --help
    <span class="go">usage: mollib measure [-h] -i id/filename [id/filename ...] [-c filename] [-l]</span>
    <span class="go">                      [-s]</span>
    <span class="go">                      [-d atom atom | -a atom atom atom | -dih atom atom atom atom | -w atom atom]</span>
    <span class="go">                      [--stats] [--only-intra] [--exclude-intra]</span>
    <span class="go">                      [--only-intra-chain] [--exclude-intra-chain]</span>
    <span class="go">                      [--only-delta DELTA] [--only-bonded] [--hydrogenate]</span>
    <span class="go">                      [--rama]</span>
    
    <span class="go">arguments:</span>
    <span class="go">  -h, --help            show this help message and exit</span>
    <span class="go">  -i id/filename [id/filename ...], --in id/filename [id/filename ...]</span>
    <span class="go">                        (required) The filename(s) or PDB identifier(s) of the</span>
    <span class="go">                        structure(s)</span>
    <span class="go">  -c filename, --config filename</span>
    <span class="go">                        The configuration filename</span>
    <span class="go">  -l                    List details on the molecule(s)</span>
    <span class="go">  -s, --save            Save fetched files to the local directory.</span>
    <span class="go">  --hydrogenate         Strip hydrogens and re-add them before analysis</span>
    
    <span class="go">measurement options:</span>
    <span class="go">  -d atom atom, --dist atom atom</span>
    <span class="go">                        Measure distances between 2 atom selections. ex: 31.N</span>
    <span class="go">                        32.CA</span>
    <span class="go">  -a atom atom atom, --angle atom atom atom</span>
    <span class="go">                        Measure angles between 3 atom selections. ex: 31.N</span>
    <span class="go">                        31.CA 31.C</span>
    <span class="go">  -dih atom atom atom atom, --dihedral atom atom atom atom</span>
    <span class="go">                        Measure dihedral angles between 4 atom selections. ex:</span>
    <span class="go">                        31.N 31.CA 31.C 32.N</span>
    <span class="go">  -w atom atom, --within atom atom</span>
    <span class="go">                        Measure all distances from atom selection to within</span>
    <span class="go">                        the specified distance. ex: 31:33.N 5</span>
    <span class="go">  --stats               Report statistics on the reported measurements.</span>
    <span class="go">  --rama                Report the Ramachandran angles. Filters and options</span>
    <span class="go">                        are ignored.</span>
    
    <span class="go">filters:</span>
    <span class="go">  --only-intra          Only report measurements within a residue</span>
    <span class="go">  --exclude-intra       Exclude measurements within a residue</span>
    <span class="go">  --only-intra-chain    Only report measurements within a chain</span>
    <span class="go">  --exclude-intra-chain</span>
    <span class="go">                        Exclude measurements within a chain</span>
    <span class="go">  --only-delta DELTA    Only report residues separated by DELTA residue</span>
    <span class="go">                        numbers</span>
    <span class="go">  --only-bonded         Only report measurements from bonded atoms</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\footnotesize]
  \textcolor{darkorange}{$} ml measure -{-}help
  usage: mollib measure [-h] -i id/filename [id/filename ...] [-c filename] [-l]
                        [-s]
                        [-d atom atom | -a atom atom atom | -dih atom atom atom atom | -w atom atom]
                        [-{-}stats] [-{-}only-intra] [-{-}exclude-intra]
                        [-{-}only-intra-chain] [-{-}exclude-intra-chain]
                        [-{-}only-delta DELTA] [-{-}only-bonded] [-{-}hydrogenate]
                        [-{-}rama]
  
  arguments:
    -h, -{-}help            show this help message and exit
    -i id/filename [id/filename ...], -{-}in id/filename [id/filename ...]
                          (required) The filename(s) or PDB identifier(s) of the
                          structure(s)
    -c filename, -{-}config filename
                          The configuration filename
    -l                    List details on the molecule(s)
    -s, -{-}save            Save fetched files to the local directory.
    -{-}hydrogenate         Strip hydrogens and re-add them before analysis
  
  measurement options:
    -d atom atom, -{-}dist atom atom
                          Measure distances between 2 atom selections. ex: 31.N
                          32.CA
    -a atom atom atom, -{-}angle atom atom atom
                          Measure angles between 3 atom selections. ex: 31.N
                          31.CA 31.C
    -dih atom atom atom atom, -{-}dihedral atom atom atom atom
                          Measure dihedral angles between 4 atom selections. ex:
                          31.N 31.CA 31.C 32.N
    -w atom atom, -{-}within atom atom
                          Measure all distances from atom selection to within
                          the specified distance. ex: 31:33.N 5
    -{-}stats               Report statistics on the reported measurements.
    -{-}rama                Report the Ramachandran angles. Filters and options
                          are ignored.
  
  filters:
    -{-}only-intra          Only report measurements within a residue
    -{-}exclude-intra       Exclude measurements within a residue
    -{-}only-intra-chain    Only report measurements within a chain
    -{-}exclude-intra-chain
                          Exclude measurements within a chain
    -{-}only-delta DELTA    Only report residues separated by DELTA residue
                          numbers
    -{-}only-bonded         Only report measurements from bonded atoms
  \end{sphinxVerbatim}
 {} 

