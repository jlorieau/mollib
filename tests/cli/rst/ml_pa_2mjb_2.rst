.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span>  ml -s pa -i 2MJB -a 2MJB --set <span class="m">0</span> --exclude CE-HE CD-HD CE-SD --fix-outliers <span class="se">\</span>
    <span class="go">  --project-methyls --summary</span>
    <span class="go"><font style="font-weight:bold;">Table: </font>Summary SVD Statistics for Molecule 2MJB</span>
    
    <span class="go">---------- --------------- ---------------- ----------------- --------------- -----------</span>
    <span class="go">Overall    Q (%): 18.7     RMS: 2.71        count: 461                                   </span>
    <span class="go">N-H        Q (%): 5.4      RMS: 0.46        count: 63         Da (Hz): 9.3    Rh: 0.147  </span>
    <span class="go">C-CA       Q (%): 17.2     RMS: 0.25        count: 58         Da (Hz): -1.6   Rh: 0.147  </span>
    <span class="go">C-H+1      Q (%): 13.8     RMS: 0.4         count: 61                                    </span>
    <span class="go">C-N+1      Q (%): 19.6     RMS: 0.17        count: 60         Da (Hz): 1.0    Rh: 0.147  </span>
    <span class="go">CA-CB      Q (%): 14.4     RMS: 0.43        count: 38         Da (Hz): -1.6   Rh: 0.147  </span>
    <span class="go">CA-HA      Q (%): 10.9     RMS: 1.9         count: 66         Da (Hz): -19.3  Rh: 0.147  </span>
    <span class="go">CB-CG      Q (%): 16.1     RMS: 2.87        count: 19                                    </span>
    <span class="go">CB-HB      Q (%): 18.3     RMS: 3.45        count: 50                                    </span>
    <span class="go">CD-CG      Q (%): 31.8     RMS: 5.67        count: 19                                    </span>
    <span class="go">CG-HG      Q (%): 43.5     RMS: 8.29        count: 27                                    </span>
    <span class="go">Alignment  Aa: 0.0004318   Ar: 6.355e-05                                                 </span>
    <span class="go">Saupe      Szz: 0.0008637  Syy: -0.0003365  Sxx: -0.0005272                              </span>
    <span class="go">Angles     Z (deg): 48.1   Y&#39; (deg): 21.8   Z&#39;&#39; (deg): 104.0                             </span>
    <span class="go">---------- --------------- ---------------- ----------------- --------------- -----------</span>
    
    <span class="go">* Inverting the sign of &#39;N-H&#39; interactions improved the overall Q-factor from 367.2% to 83.2%.</span>
    <span class="go">* Inverting the sign of &#39;C-N+1&#39; interactions improved the overall Q-factor from 83.2% to 21.1%.</span>
    <span class="go">* Removing outlier data points A.46CA-CB, A.13CB-CG2, A.16CB-HB#, A.60CB-HB#, A.24C-25N improved the</span>
    <span class="go">  overall Q-factor from 21.1% to 18.7%.</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\footnotesize]
  \textcolor{darkorange}{$}  ml -s pa -i 2MJB -a 2MJB -{-}set 0 -{-}exclude CE-HE CD-HD CE-SD -{-}fix-outliers \textbackslash
  -{-}project-methyls -{-}summary
  \textbf{Table: }Summary SVD Statistics for Molecule 2MJB
  
  -{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-
  Overall    Q (%): 18.7     RMS: 2.71        count: 461                                   
  N-H        Q (%): 5.4      RMS: 0.46        count: 63         Da (Hz): 9.3    Rh: 0.147  
  C-CA       Q (%): 17.2     RMS: 0.25        count: 58         Da (Hz): -1.6   Rh: 0.147  
  C-H+1      Q (%): 13.8     RMS: 0.4         count: 61                                    
  C-N+1      Q (%): 19.6     RMS: 0.17        count: 60         Da (Hz): 1.0    Rh: 0.147  
  CA-CB      Q (%): 14.4     RMS: 0.43        count: 38         Da (Hz): -1.6   Rh: 0.147  
  CA-HA      Q (%): 10.9     RMS: 1.9         count: 66         Da (Hz): -19.3  Rh: 0.147  
  CB-CG      Q (%): 16.1     RMS: 2.87        count: 19                                    
  CB-HB      Q (%): 18.3     RMS: 3.45        count: 50                                    
  CD-CG      Q (%): 31.8     RMS: 5.67        count: 19                                    
  CG-HG      Q (%): 43.5     RMS: 8.29        count: 27                                    
  Alignment  Aa: 0.0004318   Ar: 6.355e-05                                                 
  Saupe      Szz: 0.0008637  Syy: -0.0003365  Sxx: -0.0005272                              
  Angles     Z (deg): 48.1   Y' (deg): 21.8   Z'' (deg): 104.0                             
  -{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-
  
  * Inverting the sign of 'N-H' interactions improved the overall Q-factor from 367.2% to 83.2%.
  * Inverting the sign of 'C-N+1' interactions improved the overall Q-factor from 83.2% to 21.1%.
  * Removing outlier data points A.46CA-CB, A.13CB-CG2, A.16CB-HB#, A.60CB-HB#, A.24C-25N improved the
    overall Q-factor from 21.1% to 18.7%.
  \end{sphinxVerbatim}
 {} 

