# Differential analysis code base for Keck/HIRES observations of stars in M15 and M92

(Kirby et al. 2025, Nature, submitted)

This code base contains both IDL and Python codes.  The intial spectral preparation, the measurement of equivalent widths, and the abundance measurements per line are done in IDL.  The differential analysis, plotting, and generation of tables is done in a Python notebook.

### make_allframes.pro:

Create a binary FITS table containing paths to spectra, Gaia photometry, and atmospheric parameters derived from photometry.

### merge_all.pro:

Continuum-normalize the spectra, merge orders, and merge multiple exposures of the same star.

### hiresspec3.pro:

Graphical user interface to measure equivalent widths.

### abund.pro:

Calculate abundances for each line in a star.  Can be run in normal mode or Monte Carlo mode.

### abundall:

A bash script to run abund.pro or abundmc.pro on all stars and on different processor threads.

### runabund:

Helper function for abundall.

### m15_m92_atlas.py:

Wrapper to compute ATLAS9 model atmospheres with BasicATLAS.

### m15_m92_hires.ipynb:

A Python notebook to compute the differential abundances, then produce plots and tables.