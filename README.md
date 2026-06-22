# Enhancing-stability-in-next-generation-reservoir-computing-by-enforcing-the-echo-state-property
Here we present the matlab files that compute and plot the graphs contained in our paper

This is the code for the results and figures in our paper "Enhancing stability in next-generation reservoir computing by enforcing the echo-state property". 
They are written in Matlab, and require versions of MATLAB R2024b or R2025b.

The directories:

/Fig2/ -> plotting statistics of valid prediction time Tvp (a) and escape time from chaotic attractor Tesc (b) vs transversal Lyapunov Exponent (TLE) for Lorenz system; these results are computed for both cases: for hybrid NGRC-Ch algorithm, and for standard algorithm (SA); 

/Fig3/ -> plotting PDS and PDF for Lorenz system, computed from original Lorenz system and from hybrid NGRC-Ch system that was learned from Lorenz system;

/Fig4/ -> plotting statistics of valid prediction time Tvp and escape time from chaotic attractor Tesc vs transversal Lyapunov Exponent (TLE) for Hindmarsch Rose system; these results are computed for both cases: for hybrid NGRC-Ch algorithm, and for standard algorithm (SA);

/Fig5/ -> (a) plotting the return map interspike intervals of the membrane potential; (b) plotting the PDF of interspike intervals; 

/Fig6/ -> plotting the bifurcation diagram for original HR system (a), and for hybrid NGRC-Ch algorithm that was trained from the HR system (b); 


The auxiliary directories:

/Fig2_aux/ -> 

/Fig2_aux_Tvp/ -> Computing the statistics for Fig2(a);

/Fig2_aux_Tesc/ -> Computing statistics for Fig2(b);
		
/Fig4_aux/ ->

/Fig4_aux_Tvp/ -> Computing the statistics for Fig4(a);

/Fig4_aux_Tesc/ -> Computing statistics for Fig4(b)

/Fig6_aux/ ->

/Fig6_aux_orig/ -> Computing bifurcation diagram for Fig6(a);

/Fig6_aux_NGRC/ -> Computing bifurcation diagram for Fig6(b);

%--NOTE!!!-----------------------------------------------------------------

The figure Fig1 was plotted using the program Inkskape; the 
corresponding files are not present here;

%--------------------------------------------------------------------------

The time series (necessary for Fig3) are computed in directory 
/Fig2_aux/Fig2_aux_Tesc/;

%--------------------------------------------------------------------------

The time series (necessary for Fig5) are computed in directory 
/Fig4_aux/Fig4_aux_Tesc/;

%--------------------------------------------------------------------------

The bifurcation diagram for Fig6(b) is computed using the extended version 
of the hybrid NGRC-Ch algorithm which is described in detail in the file
[Additional_Supplementary_Material.pdf];

%--------------------------------------------------------------------------

In the directory /Supplementary_Material/ we provide the matlab files for plotting
the graphs in the Supplementary Material of our paper;


%--------------------------------------------------------------------------
