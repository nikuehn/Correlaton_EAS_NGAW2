This repository contains data and code for the paper "Bayesian interfrequency and interintensity empirical correlations of horizontal Fourier amplitude spectra using NGA-West2 data", accepted for publication in Earthquake Spectra.

'data':
  - 'TotResidAllPer_IM_ModES_combined.csv' contains total residuals of PSA and other intensity measures (PGA, CAV, AI) with respect to the C14 GMM
  - 'TotResidAllPer_EAS_ModES_combined.csv' contains total residuals of EAS and other intensity measures (PGA, CAV, AI) with respect to the C25 GMM
  - 'cor_ba19.csv' contains correlations from the BA19 analysis between a taret of 5Hz and the EAS frequencies used here.

'stan': folder contains stan models to partition total residuals for two target variables into event terms, site terms, and single-site residuals, and calculate their correlations.

'r': contains example r-scrpts to perform the analysis using cmdstanr.
 
