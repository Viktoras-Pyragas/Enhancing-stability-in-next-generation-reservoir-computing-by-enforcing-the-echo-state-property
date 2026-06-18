Comments about the MATLAB files in this directory
=================================================

This directory contains MATLAB programs for learning and testing a
next-generation reservoir computing (NGRC) model from a Kuramoto-Sivashinsky
time series. The model uses delay coordinates and Chebyshev polynomial
features to approximate the time derivative of the first spectral component
of the system.


KS_spect.m
----------
Main script of the project.

This program loads the time series from ks_spectral_a1.mat, takes the first
spectral component a(:,1), normalizes it to approximately the interval
[-1, 1], and builds a delay-coordinate representation of the signal.

The delayed state has the form

    [x(t), x(t-tau), x(t-2*tau), ..., x(t-Ndel*tau)]

The script then generates Chebyshev polynomial features from these delayed
coordinates and uses ridge regression to learn a linear model for the
derivative dx/dt.

The learned model has the general form

    dx/dt = W * r

where W is the learned coefficient vector and r is the Chebyshev feature
vector.

The script can also include a conditional Lyapunov exponent constraint using
the parameter Lc_tot. This is done by augmenting the regression system with
the derivative feature vector DR and a Lagrange multiplier.

After training, the script integrates the learned NGRC delay differential
equation forward in time, compares the predicted signal with the original
Kuramoto-Sivashinsky data, computes the valid prediction time and escape
time, produces several plots, and saves the results to files such as
Tesc_Lc.mat and XPt_record2.mat.


plot_XPt.m
----------
Plotting and post-processing script.

This program loads the saved prediction file XPt_record2.mat. The file
contains a matrix Xtot with columns corresponding to time, the exact
normalized signal, and the NGRC-predicted signal.

The script computes the normalized root-mean-square prediction error between
the exact signal and the NGRC signal. It also loads ks_spectral_a1.mat so
that it can compare the predicted dynamics with the original exact dynamics.

The main plots produced by this script are:

    1. Delay-coordinate plots of the exact system and NGRC prediction.
    2. A time-series comparison of the exact and predicted signals.

This file is useful for visual inspection of the quality of the trained NGRC
model after KS_spect.m has saved the prediction results.


generateChebyshevs.m
--------------------
Function for generating Chebyshev polynomial features from a full data set.

Input:

    X       - matrix of delayed coordinates, with one row per time point
    degree  - total polynomial degree to generate

The function first generates all exponent combinations whose sum is equal to
the chosen degree. It then computes Chebyshev polynomials of each coordinate
and forms products of these polynomials according to the exponent
combinations.

The output is a matrix where each column is one Chebyshev feature and each
row corresponds to one time point.

This function is used in KS_spect.m during training to build the feature
matrix R.


generateChebyshevsDiffLyap.m
----------------------------
Function for generating derivative-related Chebyshev features used in the
conditional Lyapunov exponent constraint.

This function is similar to generateChebyshevs.m, but it replaces the
Chebyshev polynomials of the first coordinate with their derivatives with
respect to that coordinate. It then averages the resulting feature values
over all time points.

The output DR is a column vector. In KS_spect.m, this vector is used to impose
the selected conditional Lyapunov value Lc_tot during the regression step.


generateCheb_feat_vect.m
------------------------
Function for generating a Chebyshev feature vector for one delayed state.

Input:

    X          - one row vector containing the current and delayed states
    exponents  - exponent matrix for a selected total degree
    degree     - selected total polynomial degree

The function computes Chebyshev polynomials for the single state vector and
then forms the corresponding feature products.

This function is used during the prediction/integration stage in KS_spect.m.
At each integration step, the current and delayed NGRC states are converted
into a feature vector r, and the learned model W*r gives the predicted
derivative.


generate_all_exp_matreces.m
---------------------------
Utility function for precomputing exponent matrices.

Input:

    m       - embedding dimension, equal to number of delayed coordinates
    DegMax  - maximum polynomial degree

For every degree from 1 to DegMax, the function calls genExponents.m and
stores the resulting exponent matrix in a structure called matr.

The structure is saved to matr.mat. KS_spect.m later loads this file during
the prediction stage so that exponent combinations do not need to be
regenerated at every time step.


genExponents.m
--------------
Recursive helper function for generating exponent combinations.

Input:

    m       - number of variables
    degree  - total polynomial degree

The function returns all nonnegative integer combinations of length m whose
sum is equal to degree.

For example, for m = 3 and degree = 2, the function returns combinations
such as

    [0 0 2]
    [0 1 1]
    [0 2 0]
    [1 0 1]
    [1 1 0]
    [2 0 0]

These exponent combinations define which Chebyshev polynomial order is used
for each delayed coordinate when building feature products.


Overall program flow
--------------------
The typical workflow is:

    1. Run KS_spect.m.
    2. The script loads ks_spectral_a1.mat.
    3. It builds delay coordinates from the first spectral component.
    4. It generates Chebyshev features using generateChebyshevs.m.
    5. It learns the NGRC coefficients using ridge regression.
    6. It generates exponent matrices using generate_all_exp_matreces.m.
    7. It integrates the learned NGRC model using generateCheb_feat_vect.m.
    8. It saves prediction and escape-time results.
    9. Run plot_XPt.m to inspect the saved prediction visually.

