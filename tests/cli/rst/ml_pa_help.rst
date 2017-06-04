.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span>  ml pa --help
    <span class="go">usage: mollib pa [-h] -i id/filename [id/filename ...] [-c filename] [-l] [-s]</span>
    <span class="go">                 [--hydrogenate] -a id/filename [id/filename ...]</span>
    <span class="go">                 [-o filename] [-p filename] [--summary] [--set id]</span>
    <span class="go">                 [--exclude [interaction-type [interaction-type ...]]]</span>
    <span class="go">                 [--project-methyls] [--methyl-scale number]</span>
    <span class="go">                 [--fix-sign | --nofix-sign]</span>
    <span class="go">                 [--fix-nh-scale | --nofix-nh-scale]</span>
    <span class="go">                 [--fix-outliers | --nofix-outliers]</span>
    
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
    <span class="go">  -a id/filename [id/filename ...], --alignment id/filename [id/filename ...]</span>
    <span class="go">                        (required) Alignment file or identifier with RDC and</span>
    <span class="go">                        RACS data</span>
    <span class="go">  -o filename, --out filename</span>
    <span class="go">                        The output filename for the reports of the fit data.</span>
    <span class="go">  -p filename, --pred filename</span>
    <span class="go">                        The output filename for the report of the back-</span>
    <span class="go">                        calculated RDCs and RACSs that are not in the</span>
    <span class="go">                        experimental data.</span>
    <span class="go">  --summary             Only display the fit summary</span>
    <span class="go">  --set id              If multiple datasets are available, this option</span>
    <span class="go">                        specifies which dataset to use.</span>
    <span class="go">  --exclude [interaction-type [interaction-type ...]]</span>
    <span class="go">                        Exclude one or more interactions of the following</span>
    <span class="go">                        type(s). ex: N-H or CE-HE</span>
    <span class="go">  --project-methyls     Fit methyl RDCs by projecting their values on the</span>
    <span class="go">                        corresponding C-C bond, as used by Xplor-NIH</span>
    <span class="go">  --methyl-scale number</span>
    <span class="go">                        The order parameter to use in scaling the methyl RDCs.</span>
    
    <span class="go">fixer arguments:</span>
    <span class="go">  --fix-sign            Check and fix mistakes in RDC and RACS sign</span>
    <span class="go">  --nofix-sign          Disable check in RDC and RACS sign</span>
    <span class="go">  --fix-nh-scale        Check and rescale couplings that were scaled to the</span>
    <span class="go">                        N-H RDC.</span>
    <span class="go">  --nofix-nh-scale      Disable N-H rescaling of couplings.</span>
    <span class="go">  --fix-outliers        Fit without outliers</span>
    <span class="go">  --nofix-outliers      Disable fitting without outliers</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\footnotesize]
  \textcolor{darkorange}{$}  ml pa -{-}help
  usage: mollib pa [-h] -i id/filename [id/filename ...] [-c filename] [-l] [-s]
                   [-{-}hydrogenate] -a id/filename [id/filename ...]
                   [-o filename] [-p filename] [-{-}summary] [-{-}set id]
                   [-{-}exclude [interaction-type [interaction-type ...]]]
                   [-{-}project-methyls] [-{-}methyl-scale number]
                   [-{-}fix-sign | -{-}nofix-sign]
                   [-{-}fix-nh-scale | -{-}nofix-nh-scale]
                   [-{-}fix-outliers | -{-}nofix-outliers]
  
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
    -a id/filename [id/filename ...], -{-}alignment id/filename [id/filename ...]
                          (required) Alignment file or identifier with RDC and
                          RACS data
    -o filename, -{-}out filename
                          The output filename for the reports of the fit data.
    -p filename, -{-}pred filename
                          The output filename for the report of the back-
                          calculated RDCs and RACSs that are not in the
                          experimental data.
    -{-}summary             Only display the fit summary
    -{-}set id              If multiple datasets are available, this option
                          specifies which dataset to use.
    -{-}exclude [interaction-type [interaction-type ...]]
                          Exclude one or more interactions of the following
                          type(s). ex: N-H or CE-HE
    -{-}project-methyls     Fit methyl RDCs by projecting their values on the
                          corresponding C-C bond, as used by Xplor-NIH
    -{-}methyl-scale number
                          The order parameter to use in scaling the methyl RDCs.
  
  fixer arguments:
    -{-}fix-sign            Check and fix mistakes in RDC and RACS sign
    -{-}nofix-sign          Disable check in RDC and RACS sign
    -{-}fix-nh-scale        Check and rescale couplings that were scaled to the
                          N-H RDC.
    -{-}nofix-nh-scale      Disable N-H rescaling of couplings.
    -{-}fix-outliers        Fit without outliers
    -{-}nofix-outliers      Disable fitting without outliers
  \end{sphinxVerbatim}
 {} 

