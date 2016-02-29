# STA-LTA-OS: STA-LTA with outlier statistics

## Introduction
This is the online repository for code from `[1]`. Code is stable, semi-flexible, and has been cleaned and optimized from the original scripts.

# Aim
STA-LTA with outlier statistics improves detection and picking accuracy, particularly in environments with many closely spaced events (e.g. earthquake swarms, hydraulic fractures). 

# Approach
In these programs, STA-LTA is recast as an outlier statistics problem using a two-state hidden Markov model (HMM), which I resolve with an expectation-maximization (EM) algorithm to determine the state membership probabilities of transformed data points.

# Use
This program can be used without modification by creating a matrix **X** of column vector data and typing the command 
`Ev = hmmscan(X);`

# Descriptions of Individual Programs
* **hmmscan**: Scan a record of data by applying EM to successive long data windows.
* **exmax**: Two-state classification of column vector data using an expectation-maximization algorithm. 
* **genprep**: Generic data preprocessing.
* **es94ap**: An automatic P-picker based on `[2]`. 
* **hmminit**: Set initial HMM parameters.
* **sop**: Event detection and picking from short-term averages of outlier probabilities.

# Limitations
1. Because **hmmscan** expects a matrix of data in column vectors, it's implicitly assumed that data are sampled at a uniform rate.
2. The default program options are intended for borehole data sampled at 4 kHz.

# Supplemental Notes
1. Program **exmax** contains two solvers that weren't used in the paper: An ordinary Gaussian solver (opts.slv = 'gauss') and a questionable (i.e. possibly wrong) solver for beta distributions (opts.slv = 'beta'). I only tested the beta solver to the extent that I can confidently say it usually converges when fed a time series of beta-distributed coefficients. Use at your peril.
2. I do not have data distribution rights. Data from the paper can be obtained from professor Mirko van der Baan (University of Alberta) or professor David Eaton (University of Calgary).

## References
1. Jones, J.P., and van der Baan, M. (2015). Adaptive STA-LTA with outlier statistics, Bull. Seismol. Soc. Am.  105 (3), 1606--1618. doi: 10.1785/0120140203.
2. Earle, P., and P. Shearer (1994). Characterization of global seismograms using an automatic picking algorithm, Bull. Seismol. Soc. Am. 84 (2), 366--376.
