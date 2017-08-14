.. only:: html

.. raw:: html

    <div class="highlight"><pre><span></span><span class="gp">$</span>  ml -s pa -i 1UBQ -d 2MJB --set <span class="m">0</span> --fix-outliers --project-methyls <span class="se">\</span>
    <span class="go">  --hydrogenate --summary</span>
    <span class="go"><font style="font-weight:bold;">Table: </font>Summary SVD Statistics for Molecule 1UBQ</span>
    
    <span class="go">---------- --------------- --------------- ----------------- --------------- -----------</span>
    <span class="go">Overall    Q (%): 39.8     RMS: 5.77       count: 474                                   </span>
    <span class="go">N-H        Q (%): 14.7     RMS: 1.2        count: 62         Da (Hz): 8.9    Rh: 0.181  </span>
    <span class="go">C-CA       Q (%): 26.2     RMS: 0.37       count: 58         Da (Hz): -1.5   Rh: 0.181  </span>
    <span class="go">C-H+1      Q (%): 20.5     RMS: 0.56       count: 60                                    </span>
    <span class="go">C-N+1      Q (%): 31.1     RMS: 0.26       count: 61         Da (Hz): 0.9    Rh: 0.181  </span>
    <span class="go">CA-CB      Q (%): 25.1     RMS: 0.36       count: 37         Da (Hz): -1.5   Rh: 0.181  </span>
    <span class="go">CA-HA      Q (%): 25.8     RMS: 4.32       count: 66         Da (Hz): -18.4  Rh: 0.181  </span>
    <span class="go">CB-CG      Q (%): 24.4     RMS: 4.18       count: 18                                    </span>
    <span class="go">CB-HB      Q (%): 63.8     RMS: 10.42      count: 51                                    </span>
    <span class="go">CD-CG      Q (%): 55.1     RMS: 9.42       count: 19                                    </span>
    <span class="go">CD-HD      Q (%): 49.6     RMS: 8.45       count: 10                                    </span>
    <span class="go">CE-HE      Q (%): 19.7     RMS: 3.68       count: 4                                     </span>
    <span class="go">CE-SD      Q (%): 36.7     RMS: -          count: 1                                     </span>
    <span class="go">CG-HG      Q (%): 94.9     RMS: 15.63      count: 27                                    </span>
    <span class="go">Alignment  Aa: 0.000412    Ar: 7.466e-05                                                </span>
    <span class="go">Saupe      Szz: 0.0008239  Syy: -0.0003    Sxx: -0.0005239                              </span>
    <span class="go">Angles     Z (deg): 330.1  Y&#39; (deg): 70.4  Z&#39;&#39; (deg): 257.7                             </span>
    <span class="go">---------- --------------- --------------- ----------------- --------------- -----------</span>
    
    <span class="go">* Inverting the sign of &#39;N-H&#39; interactions improved the overall Q-factor from 516.4% to 93.7%.</span>
    <span class="go">* Inverting the sign of &#39;C-N+1&#39; interactions improved the overall Q-factor from 93.7% to 46.3%.</span>
    <span class="go">* Removing outlier data points A.46CA-CB, A.14CB-CG2, A.48CD-HD#, A.33CE-HE#, A.48N-H, A.7C-8H,</span>
    <span class="go">  A.28CA-CB, A.13CB-CG2, A.14CB-HB improved the overall Q-factor from 46.3% to 39.8%.</span>
    </pre></div>


.. only:: latex

.. raw:: latex

  \begin{sphinxVerbatim}[commandchars=\\\{\},fontsize=\small]
  \textcolor{darkorange}{$}  ml -s pa -i 1UBQ -d 2MJB -{-}set 0 -{-}fix-outliers -{-}project-methyls \textbackslash
  -{-}hydrogenate -{-}summary
  \textbf{Table: }Summary SVD Statistics for Molecule 1UBQ
  
  -{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-
  Overall    Q (%): 39.8     RMS: 5.77       count: 474                                   
  N-H        Q (%): 14.7     RMS: 1.2        count: 62         Da (Hz): 8.9    Rh: 0.181  
  C-CA       Q (%): 26.2     RMS: 0.37       count: 58         Da (Hz): -1.5   Rh: 0.181  
  C-H+1      Q (%): 20.5     RMS: 0.56       count: 60                                    
  C-N+1      Q (%): 31.1     RMS: 0.26       count: 61         Da (Hz): 0.9    Rh: 0.181  
  CA-CB      Q (%): 25.1     RMS: 0.36       count: 37         Da (Hz): -1.5   Rh: 0.181  
  CA-HA      Q (%): 25.8     RMS: 4.32       count: 66         Da (Hz): -18.4  Rh: 0.181  
  CB-CG      Q (%): 24.4     RMS: 4.18       count: 18                                    
  CB-HB      Q (%): 63.8     RMS: 10.42      count: 51                                    
  CD-CG      Q (%): 55.1     RMS: 9.42       count: 19                                    
  CD-HD      Q (%): 49.6     RMS: 8.45       count: 10                                    
  CE-HE      Q (%): 19.7     RMS: 3.68       count: 4                                     
  CE-SD      Q (%): 36.7     RMS: -          count: 1                                     
  CG-HG      Q (%): 94.9     RMS: 15.63      count: 27                                    
  Alignment  Aa: 0.000412    Ar: 7.466e-05                                                
  Saupe      Szz: 0.0008239  Syy: -0.0003    Sxx: -0.0005239                              
  Angles     Z (deg): 330.1  Y' (deg): 70.4  Z'' (deg): 257.7                             
  -{-}-{-}-{-}-{-}-{-} -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-{-}-{-}- -{-}-{-}-{-}-{-}-{-}-
  
  * Inverting the sign of 'N-H' interactions improved the overall Q-factor from 516.4% to 93.7%.
  * Inverting the sign of 'C-N+1' interactions improved the overall Q-factor from 93.7% to 46.3%.
  * Removing outlier data points A.46CA-CB, A.14CB-CG2, A.48CD-HD#, A.33CE-HE#, A.48N-H, A.7C-8H,
    A.28CA-CB, A.13CB-CG2, A.14CB-HB improved the overall Q-factor from 46.3% to 39.8%.
  \end{sphinxVerbatim}
 {} 

