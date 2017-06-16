.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span>  ml process --help
    <span class="go">usage: mollib process [-h] -i id/filename [id/filename ...] [-c filename] [-l]</span>
    <span class="go">                      [-s] [-m [MODELS [MODELS ...]]]</span>
    <span class="go">                      [-o [filename [filename ...]]] [--hydrogenate]</span>
    
    <span class="go">arguments:</span>
    <span class="go">  -h, --help            show this help message and exit</span>
    <span class="go">  -i id/filename [id/filename ...], --in id/filename [id/filename ...]</span>
    <span class="go">                        (required) The filename(s) or PDB identifier(s) of the</span>
    <span class="go">                        structure(s)</span>
    <span class="go">  -c filename, --config filename</span>
    <span class="go">                        The configuration filename</span>
    <span class="go">  -l                    List details on the molecule(s)</span>
    <span class="go">  -s, --save            Save fetched files to the local directory.</span>
    <span class="go">  -m [MODELS [MODELS ...]], --models [MODELS [MODELS ...]]</span>
    <span class="go">                        The models numbers to analyze.</span>
    <span class="go">  -o [filename [filename ...]], --out [filename [filename ...]]</span>
    <span class="go">                        The output filename(s) for the structure(s)</span>
    <span class="go">  --hydrogenate         Strip hydrogens and re-add them before analysis</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\small]
  \textcolor{darkorange}{$}  ml process -{-}help
  usage: mollib process [-h] -i id/filename [id/filename ...] [-c filename] [-l]
                        [-s] [-m [MODELS [MODELS ...]]]
                        [-o [filename [filename ...]]] [-{-}hydrogenate]
  
  arguments:
    -h, -{-}help            show this help message and exit
    -i id/filename [id/filename ...], -{-}in id/filename [id/filename ...]
                          (required) The filename(s) or PDB identifier(s) of the
                          structure(s)
    -c filename, -{-}config filename
                          The configuration filename
    -l                    List details on the molecule(s)
    -s, -{-}save            Save fetched files to the local directory.
    -m [MODELS [MODELS ...]], -{-}models [MODELS [MODELS ...]]
                          The models numbers to analyze.
    -o [filename [filename ...]], -{-}out [filename [filename ...]]
                          The output filename(s) for the structure(s)
    -{-}hydrogenate         Strip hydrogens and re-add them before analysis
  \end{sphinxVerbatim}
 {} 

