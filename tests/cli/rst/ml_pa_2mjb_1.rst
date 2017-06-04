.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span>  ml pa -i 2MJB -a 2MJB --set <span class="m">0</span> --fix-outliers --project-methyls --summary
    <span class="go"><font style="font-weight:bold;">Table: </font>Summary SVD Statistics for Molecule 2MJB</span>
    
    <span class="go">---------- --------------- ---------------- ----------------- --------------- -----------</span>
    <span class="go">Overall    Q (%): 23.8     RMS: 3.68        count: 477                                   </span>
    <span class="go">N-H        Q (%): 6.3      RMS: 0.52        count: 63         Da (Hz): 9.1    Rh: 0.144  </span>
    <span class="go">C-CA       Q (%): 19.4     RMS: 0.28        count: 58         Da (Hz): -1.6   Rh: 0.144  </span>
    <span class="go">C-H+1      Q (%): 13.0     RMS: 0.37        count: 61                                    </span>
    <span class="go">C-N+1      Q (%): 22.1     RMS: 0.19        count: 60         Da (Hz): 0.9    Rh: 0.144  </span>
    <span class="go">CA-CB      Q (%): 15.0     RMS: 0.4         count: 38         Da (Hz): -1.6   Rh: 0.144  </span>
    <span class="go">CA-HA      Q (%): 12.9     RMS: 2.21        count: 66         Da (Hz): -18.8  Rh: 0.144  </span>
    <span class="go">CB-CG      Q (%): 16.8     RMS: 2.93        count: 19                                    </span>
    <span class="go">CB-HB      Q (%): 18.1     RMS: 3.34        count: 50                                    </span>
    <span class="go">CD-CG      Q (%): 31.2     RMS: 5.44        count: 19                                    </span>
    <span class="go">CD-HD      Q (%): 49.5     RMS: 9.53        count: 10                                    </span>
    <span class="go">CE-HE      Q (%): 115.0    RMS: 23.44       count: 5                                     </span>
    <span class="go">CE-SD      Q (%): 74.1     RMS: -           count: 1                                     </span>
    <span class="go">CG-HG      Q (%): 43.1     RMS: 8.04        count: 27                                    </span>
    <span class="go">Alignment  Aa: 0.0004225   Ar: 6.07e-05                                                  </span>
    <span class="go">Saupe      Szz: 0.0008451  Syy: -0.0003315  Sxx: -0.0005136                              </span>
    <span class="go">Angles     Z (deg): 49.0   Y&#39; (deg): 22.2   Z&#39;&#39; (deg): 104.8                             </span>
    <span class="go">---------- --------------- ---------------- ----------------- --------------- -----------</span>
    
    <span class="go">* Inverting the sign of &#39;N-H&#39; interactions improved the overall Q-factor from 351.1% to 84.6%.</span>
    <span class="go">* Inverting the sign of &#39;C-N+1&#39; interactions improved the overall Q-factor from 84.6% to 27.1%.</span>
    <span class="go">* Removing outlier data points A.46CA-CB, A.48CD-HD#, A.16CB-HB#, A.60CB-HB#, A.13CB-CG2, A.24C-25N</span>
    <span class="go">  improved the overall Q-factor from 27.1% to 23.8%.</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\footnotesize]
  \textcolor{darkorange}{$}  ml pa -i 2MJB -a 2MJB -{-}set 0 -{-}fix-outliers -{-}project-methyls -{-}summary
  \textbf{Table: }Summary SVD Statistics for Molecule 2MJB
  
  -{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-
  Overall    Q (%): 23.8     RMS: 3.68        count: 477                                   
  N-H        Q (%): 6.3      RMS: 0.52        count: 63         Da (Hz): 9.1    Rh: 0.144  
  C-CA       Q (%): 19.4     RMS: 0.28        count: 58         Da (Hz): -1.6   Rh: 0.144  
  C-H+1      Q (%): 13.0     RMS: 0.37        count: 61                                    
  C-N+1      Q (%): 22.1     RMS: 0.19        count: 60         Da (Hz): 0.9    Rh: 0.144  
  CA-CB      Q (%): 15.0     RMS: 0.4         count: 38         Da (Hz): -1.6   Rh: 0.144  
  CA-HA      Q (%): 12.9     RMS: 2.21        count: 66         Da (Hz): -18.8  Rh: 0.144  
  CB-CG      Q (%): 16.8     RMS: 2.93        count: 19                                    
  CB-HB      Q (%): 18.1     RMS: 3.34        count: 50                                    
  CD-CG      Q (%): 31.2     RMS: 5.44        count: 19                                    
  CD-HD      Q (%): 49.5     RMS: 9.53        count: 10                                    
  CE-HE      Q (%): 115.0    RMS: 23.44       count: 5                                     
  CE-SD      Q (%): 74.1     RMS: -           count: 1                                     
  CG-HG      Q (%): 43.1     RMS: 8.04        count: 27                                    
  Alignment  Aa: 0.0004225   Ar: 6.07e-05                                                  
  Saupe      Szz: 0.0008451  Syy: -0.0003315  Sxx: -0.0005136                              
  Angles     Z (deg): 49.0   Y' (deg): 22.2   Z'' (deg): 104.8                             
  -{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-
  
  * Inverting the sign of 'N-H' interactions improved the overall Q-factor from 351.1% to 84.6%.
  * Inverting the sign of 'C-N+1' interactions improved the overall Q-factor from 84.6% to 27.1%.
  * Removing outlier data points A.46CA-CB, A.48CD-HD#, A.16CB-HB#, A.60CB-HB#, A.13CB-CG2, A.24C-25N
    improved the overall Q-factor from 27.1% to 23.8%.
  \end{sphinxVerbatim}
 {} 

