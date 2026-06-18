Comments about MATLAB m-files in this directory

MG_Lc_aut.m
Main script for learning and testing an NGRC delay-differential model from
Mackey-Glass data. It loads the MG time series, normalizes the training
signal, builds a delayed embedding, constructs Chebyshev feature matrices,
solves a ridge-regression problem with a conditional Lyapunov constraint
controlled by Lc, and then integrates the learned closed NGRC system. The
script computes prediction/escape diagnostics, saves Xtot_record.mat, and
produces comparison plots between the true Mackey-Glass signal and the NGRC
replication.

generateChebyshevs.m
Function that builds multivariate Chebyshev feature columns of one specified
total degree from an input data matrix X. Each row of X is a time point and
each column is an embedded coordinate. The function generates all exponent
combinations whose sum equals the requested degree, evaluates Chebyshev
polynomials for each coordinate, and returns the corresponding monomial-like
feature matrix.

generateChebyshevsDiffLyap.m
Function that computes the averaged derivative-related Chebyshev feature
vector used by MG_Lc_aut.m to impose the conditional Lyapunov constraint. It
generates the same total-degree exponent combinations as generateChebyshevs.m,
but replaces the first coordinate's Chebyshev terms with derivatives before
forming products and averaging each feature over all time points.

generateCheb_feat_vect.m
Function for evaluating one block of the Chebyshev feature vector at a single
embedded state X. It receives a precomputed exponent matrix and a degree,
constructs one-dimensional Chebyshev polynomials for each embedded coordinate,
and returns the feature values needed during NGRC time integration.

generate_all_exp_matreces.m
Helper function that precomputes exponent matrices for all polynomial degrees
from 1 through DegMax for an embedding dimension m. It stores the matrices in a
struct array named matr and saves them to matr.mat, which MG_Lc_aut.m later
loads during the prediction loop.

genExponents.m
Recursive helper function that returns every length-m nonnegative integer
exponent combination whose entries sum to the requested total degree. This is
the shared combinatorial basis used for arranging Chebyshev feature products.

plot_MG_rc.m
Small plotting script for visualizing a delayed-coordinate phase portrait from
Xtot_record.mat. It loads the saved true/predicted time series, forms the pair
x(t) and x(t + tau) from the true signal, truncates the plot to a fixed window,
and draws the resulting delay embedding curve.

