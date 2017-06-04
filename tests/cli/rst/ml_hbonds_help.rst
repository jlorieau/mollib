.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span>  ml hbonds --help
    <span class="go">usage: mollib hbonds [-h] -i id/filename [id/filename ...] [-c filename] [-l]</span>
    <span class="go">                     [-s] [--hydrogenate] [--aliphatic] [--detailed]</span>
    <span class="go">                     [--sort-type]</span>
    
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
    
    <span class="go">hbond options:</span>
    <span class="go">  --aliphatic           Includes aliphatic hydrogen bonds</span>
    <span class="go">  --detailed            Report detailed information on hydrogen bonds.</span>
    <span class="go">  --sort-type           Sort hydrogen bonds by type</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\footnotesize]
  \textcolor{darkorange}{$}  ml hbonds -{-}help
  usage: mollib hbonds [-h] -i id/filename [id/filename ...] [-c filename] [-l]
                       [-s] [-{-}hydrogenate] [-{-}aliphatic] [-{-}detailed]
                       [-{-}sort-type]
  
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
  
  hbond options:
    -{-}aliphatic           Includes aliphatic hydrogen bonds
    -{-}detailed            Report detailed information on hydrogen bonds.
    -{-}sort-type           Sort hydrogen bonds by type
  \end{sphinxVerbatim}
 {} 

