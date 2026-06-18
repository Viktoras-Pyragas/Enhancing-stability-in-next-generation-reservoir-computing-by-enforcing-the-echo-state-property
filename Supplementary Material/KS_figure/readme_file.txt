README FOR THE MATLAB FILES IN THIS DIRECTORY
=============================================

This directory contains MATLAB code for plotting and analyzing results for
the Kuramoto-Sivashinsky system. The scripts compare an original/exact time
series with an autonomous NGRC reservoir prediction.


Main script
-----------

main_Fig_KS.m

This is the main plotting script. It loads the data file Xtot_record.mat,
extracts the time vector, the exact/original solution, and the NGRC predicted
solution, and creates the final comparison figure.

The script assumes that the loaded variable Xtot has three columns:

  Xtot(:,1) = Tpr1  - time samples
  Xtot(:,2) = U1    - original/exact system dynamics
  Xtot(:,3) = XP1   - autonomous NGRC reservoir dynamics

The script produces a multi-panel figure:

  (a) Time series comparison of the original and NGRC signals.
  (b) Delayed-coordinate phase portrait of the original signal.
  (c) Delayed-coordinate phase portrait of the NGRC signal.
  (d) Probability density function comparison.
  (e) Power spectral density comparison.

At the end, the script saves the generated figure as Fig_KS.eps.


Helper functions
----------------

plot_KS_dyn.m

This function prepares delayed-coordinate data for the phase portraits.
It receives Xtot as input and extracts the original signal U1 and the NGRC
signal XP1.

It also calculates and prints the valid prediction time, defined as the first
time when the absolute difference between the NGRC prediction and the original
signal becomes larger than 0.02.

The function uses the delay tau = 0.17 and returns:

  Xnow     - NGRC signal at time t
  Xdel     - NGRC signal at delayed time t - tau
  Xnow_ex  - original signal at time t
  Xdel_ex  - original signal at delayed time t - tau


plot_KS_pdf.m

This function computes probability density functions for the original and
NGRC signals.

It removes non-finite values such as NaN or Inf, uses common histogram bin
limits for both signals, and normalizes the histograms as probability density
functions. Using common bins makes the two PDFs directly comparable.

The function returns:

  bin_centers - centers of the histogram bins
  pdf_ngrc    - PDF of the NGRC signal
  pdf_exact   - PDF of the original/exact signal


plot_KS_psd.m

This function computes the power spectral density of the original and NGRC
signals using Welch's method.

It extracts the time vector and both signals from Xtot, removes invalid data
points, computes the sampling frequency from the time vector, subtracts the
mean from each signal, and then applies MATLAB's pwelch function.

The function returns:

  PSD_exact - power spectral density of the original/exact signal
  PSD_ngrc  - power spectral density of the NGRC signal
  f         - frequency vector

In the main script, the PSD values are currently loaded from KS_psd.mat instead
of being recomputed every time. The commented lines in the main script show how
to recompute and save the PSD data if needed.


Data files used by the MATLAB code
----------------------------------

Xtot_record.mat

This file contains the main data matrix Xtot used by the plotting scripts.
The expected column order is time, original/exact signal, and NGRC prediction.


KS_psd.mat

This file contains precomputed spectral data:

  f
  PSD_exact
  PSD_ngrc

The main script loads this file for panel (e), which avoids recalculating the
power spectral density on every run.


Typical use
-----------

To generate the figure, run the main script in MATLAB:

  main_Fig_KS

The helper functions are called automatically by the main script. They do not
need to be run separately unless you want to inspect one specific calculation.

