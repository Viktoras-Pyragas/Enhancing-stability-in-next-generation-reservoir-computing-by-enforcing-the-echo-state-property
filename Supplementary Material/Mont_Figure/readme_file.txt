MATLAB FILE COMMENTS
====================

This directory contains MATLAB programs for comparing an original signal with
an NGRC/reservoir-computing reconstructed signal. The input data are loaded
from TUV.mat, where the matrix TUV is expected to contain:

  column 1: time
  column 2: original signal U
  column 3: predicted/reconstructed signal V


mainT.m
-------------
Main driver function for the whole analysis.

This file:
  - clears the workspace and closes old figures;
  - defines plotting parameters such as colors, line widths, and time windows;
  - loads TUV.mat;
  - packages parameters using setParam.m;
  - packages data using setData.m;
  - calls Plot_Short_Long_Intersp_map.m to create time-domain, phase-space,
    and interspike-interval plots;
  - calls Plot_spectra.m to create the power spectrum comparison;
  - saves the resulting figure as Fig_Mnt.fig and Fig_Mnt.eps;
  - prints FIN when complete.

Important parameters in this file:
  Dt    = total time interval used for analysis and spectra;
  Nwind = number of windows used in Welch spectral analysis;
  TN    = short-time interval shown in the first plot;
  TL    = long-time interval used for phase portraits.


setParam.m
----------
Utility function that stores plotting and analysis parameters in one structure.

Input:
  individual parameters such as line widths, colors, transparency controls,
  segment length, analysis time interval, number of Welch windows, and plot
  durations.

Output:
  Qp, a structure containing all parameter values.

Purpose:
  This avoids passing many separate arguments into every plotting function.


getParam.m
----------
Utility function that extracts the individual parameter values from Qp.

Input:
  Qp, the parameter structure created by setParam.m.

Output:
  the original individual values:
  lw2, lw3, lwd, color_ex, color_rc, a1, a2, segLength, Dt, Nwind, TN, TL.

Purpose:
  This is the inverse helper of setParam.m.


setData.m
---------
Utility function that stores the TUV data matrix in one structure.

Input:
  TUV, the matrix containing time, original signal, and predicted signal.

Output:
  Qd, a structure with field Qd.TUV.

Purpose:
  This gives the plotting functions a simple way to receive the data.


getData.m
---------
Utility function that extracts the TUV matrix from the Qd data structure.

Input:
  Qd, the structure created by setData.m.

Output:
  TUV, the original data matrix.

Purpose:
  This is the inverse helper of setData.m.


Plot_Short_Long_Intersp_map.m
-----------------------------
Creates the main time-domain and nonlinear-dynamics plots.

This function:
  - reads parameters from Qp using getParam.m;
  - reads signal data from Qd using getData.m;
  - computes the time step h from the first two time values;
  - keeps only the first Dt milliseconds of data;
  - separates TUV into time TT, original signal U, and predicted signal V;
  - detects local maxima of U and V using islocalmax;
  - keeps positive maxima for interspike-interval analysis.

Generated subplots:

  (a) Short-term signal comparison
      Plots original signal U and predicted signal V over the first TN ms.
      This shows how closely the prediction follows the original signal in
      the time domain.

  (b) Phase portrait of the original signal
      Plots u(t) versus u(t - tau), where tau = 0.6 ms.
      This shows the attractor-like structure of the original dynamics.

  (c) Phase portrait of the predicted signal
      Plots v(t) versus v(t - tau), using the same delay as subplot (b).
      This helps compare whether the predicted signal reproduces the same
      dynamical structure.

  (d) Interspike interval map
      Computes time intervals between successive positive peaks and plots
      T_n versus T_{n+1}. Original and predicted signals are shown together.
      This checks whether spike timing dynamics are preserved.


Plot_spectra.m
--------------
Creates the power spectral density comparison plot.

This function:
  - reads parameters from Qp using getParam.m;
  - reads signal data from Qd using getData.m;
  - computes the sampling frequency from the time step;
  - creates a Hann window with length based on the number of Welch windows;
  - calls WelchPowerSpectralDensity.m for the original signal U;
  - calls WelchPowerSpectralDensity.m for the predicted signal V;
  - plots both spectra on a logarithmic power scale.

Generated subplot:

  (e) Power spectrum comparison
      Compares the frequency content of the original and predicted signals
      from 0 to 2 kHz. Similar spectra indicate that the prediction preserves
      important frequency-domain properties of the original signal.


Hann.m
------
Creates a Hann window of a specified integer length.

Input:
  Length, the desired number of points in the window.

Output:
  Hann, a column vector containing the Hann window values.

Purpose:
  The Hann window is used before FFT calculations in Welch spectral analysis.
  Windowing reduces spectral leakage when finite signal segments are analyzed.


WelchPowerSpectralDensity.m
---------------------------
Computes a power spectral density estimate using Welch's method.

Inputs:
  Data1             primary input signal;
  Data2             optional second signal for cross-spectrum calculation;
                    use [] for ordinary auto-spectrum calculation;
  Window            window vector, usually created by Hann.m;
  OverlapFactor     fraction of overlap between neighboring windows;
  SamplingFrequency sampling frequency of the signal.

Outputs:
  Power             estimated power spectrum;
  Frequency         frequency values corresponding to Power.

Algorithm:
  - checks that input data and window values are valid vectors;
  - divides the signal into windowed segments;
  - applies the window to each segment;
  - computes the FFT of each segment;
  - accumulates and averages the power values;
  - normalizes the result using the window energy and sampling frequency;
  - returns only the positive-frequency half of the spectrum.

Purpose:
  This function is the spectral-analysis engine used by Plot_spectra.m.


Overall workflow
----------------
Run:

  mainT

The program loads TUV.mat, creates the comparison figure, saves the figure
files, and displays FIN when the analysis is complete.

