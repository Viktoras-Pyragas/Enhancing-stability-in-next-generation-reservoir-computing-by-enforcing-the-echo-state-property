%-----------------------------------------------------------------------
% Comments about matlab files involved in this directory
%-----------------------------------------------------------------------

QIF_Cheb_lam.m
Main script for nonlinear prediction of the QIF network signal. It loads the
QIF data, selects one signal component, rescales it to approximately [-1, 1],
builds a delayed-coordinate embedding, and fits a Chebyshev-polynomial model
to the signal derivative. The regression includes ridge regularization and an
extra constraint that forces the reconstructed model to have a prescribed
conditional Lyapunov exponent, set by the variable lam. After training, the
script iterates the learned model forward, compares the forecast with the real
signal in a plot, and saves [time, real signal, predicted signal] to TUV.mat.

generateChebyshevs.m
Function that generates Chebyshev-polynomial feature columns for a full data
matrix X. Each row of X is one embedded state, and each column is one delayed
coordinate. For a requested total polynomial degree, the function builds all
products of Chebyshev polynomials whose orders sum to that degree. These
features are used in QIF_Cheb_lam.m during model fitting.

generateChebyshevsDiffLyap.m
Function that generates averaged derivatives of the Chebyshev feature terms
with respect to the first embedding coordinate. These derivative averages form
the vector DR used in QIF_Cheb_lam.m to impose the conditional Lyapunov
exponent constraint W * DR = lam during regression.

generateCheb_feat_vect.m
Function that generates a Chebyshev feature vector for a single embedded state
instead of a whole data matrix. It is used inside the prediction loop in
QIF_Cheb_lam.m, where the learned model must evaluate the nonlinear feature
vector at each forecast step.

generate_all_exp_matreces.m
Utility function that precomputes exponent-combination matrices for all
polynomial degrees from 1 to DegMax. The matrices are stored in a struct called
matr and saved to matr.mat. QIF_Cheb_lam.m loads this file during prediction so
that exponent combinations do not need to be regenerated repeatedly.

genExponents.m
Recursive utility function that generates all nonnegative integer exponent
combinations of length m whose entries sum to a chosen polynomial degree. These
combinations define which Chebyshev polynomial order is applied to each
embedding coordinate in a product feature.

Hamming.m
Function that returns a Hamming window of a requested integer length. This is a
standard signal-processing window, useful before FFT or spectral-density
calculations to reduce boundary artifacts.

Hann.m
Function that returns a Hann window of a requested integer length. Like the
Hamming window, it is commonly used in spectral analysis before applying an FFT.

WelchPowerSpectralDensity.m
Function that estimates a power spectral density using Welch's method. It can
compute the spectrum of one signal or the cross-spectrum of two signals. The
input data are split into overlapping, windowed segments; each segment is
transformed by FFT; and the resulting spectra are averaged. The function returns
the estimated power and the corresponding frequency vector.

