README for MATLAB files in this directory
=========================================

This directory contains MATLAB programs for comparing the original
Mackey-Glass system with an autonomous NGRC/reservoir prediction.


1. main_Fig_MG.m
----------------------

This is the main script. Run this file to generate the final multi-panel
figure.

The script loads Xtot_record.mat. The loaded matrix Xtot is expected to
contain three columns:

  Column 1: Tpr1 - time samples
  Column 2: U1   - original Mackey-Glass time series
  Column 3: XP1  - autonomous NGRC/reservoir time series

The script creates a figure with five panels:

  (a) Time-domain comparison of original and NGRC signals.
  (b) Delayed phase portrait of the original Mackey-Glass system.
  (c) Delayed phase portrait of the NGRC/reservoir system.
  (d) Probability density function comparison.
  (e) Power spectral density comparison.

The script calls the three helper functions:

  plot_MG_dyn.m
  plot_MG_pdf.m
  plot_MG_psd.m

At the end, it saves the final figure as Fig_MG.eps.


2. plot_MG_dyn.m
----------------

This function prepares delayed-coordinate data for the phase portraits.

Function call:

  [Xnow,Xdel,Xnow_ex,Xdel_ex] = plot_MG_dyn(Xtot)

Input:

  Xtot - matrix containing time, original signal, and NGRC signal.

Outputs:

  Xnow    - NGRC signal at the current time.
  Xdel    - NGRC signal delayed by tau.
  Xnow_ex - original signal at the current time.
  Xdel_ex - original signal delayed by tau.

The function also estimates the valid prediction time. It compares the
original and NGRC signals and finds the first point where the absolute
difference is larger than 0.02. This time is printed in the MATLAB command
window as Tvp.

The delay used for the phase portrait is:

  tau1 = 2.0


3. plot_MG_pdf.m
----------------

This function computes probability density functions for the original and
NGRC signals.

Function call:

  [bin_centers,pdf_ngrc,pdf_exact] = plot_MG_pdf(Xtot)

Input:

  Xtot - matrix containing time, original signal, and NGRC signal.

Outputs:

  bin_centers - centers of the histogram bins.
  pdf_ngrc    - estimated PDF of the NGRC/reservoir signal.
  pdf_exact   - estimated PDF of the original Mackey-Glass signal.

The function removes NaN and Inf values before computing the histograms. It
uses the same bin edges for both signals so that the two PDFs can be compared
directly.


4. plot_MG_psd.m
----------------

This function computes the power spectral density of the original and NGRC
signals using Welch's method.

Function call:

  [PSD_exact,PSD_ngrc,f] = plot_MG_psd(Xtot)

Input:

  Xtot - matrix containing time, original signal, and NGRC signal.

Outputs:

  PSD_exact - power spectral density of the original Mackey-Glass signal.
  PSD_ngrc  - power spectral density of the NGRC/reservoir signal.
  f         - frequency vector.

The function removes invalid values, estimates the sampling frequency from
the time vector, subtracts the mean from both signals, and then calls pwelch
to estimate the spectra.


General workflow
----------------

1. Make sure Xtot_record.mat is in the same directory.
2. Run main_Fig_MG.m in MATLAB.
3. The script will load the data, call the helper functions, display the
   figure, and save the result as Fig_MG.eps.

