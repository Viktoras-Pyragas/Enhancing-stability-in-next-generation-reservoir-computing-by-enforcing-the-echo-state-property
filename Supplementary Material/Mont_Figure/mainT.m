function mainT()
% mainT
% Main driver for the NGRC signal comparison plots.
%
% The script loads time-series data from TUV.mat, where TUV contains time,
% the original signal, and the NGRC reconstructed signal. It then creates a
% multi-panel figure comparing the signals in the time domain, phase space,
% interspike-interval map, and frequency domain.

close all
clearvars
format long
clc

oldFigureVisible = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
cleanupFigureVisible = onCleanup(@() set(groot, 'DefaultFigureVisible', oldFigureVisible));

% Plot and analysis parameters
lwd = 0.7; % line width for phase portraits, Figs. (b) and (c)
lw2 = 0.2; % line width for time-series and spectral plots, Figs. (a) and (e)
lw3 = 0.5; % marker size for the interspike-interval map, Fig. (d)
Dt = 30000; % time interval used for the interspike map and PSD calculation
Nwind=128; % number of windows used in Welch spectral estimation
TN = 800; % time interval shown in the short-term dynamics plot, Fig. (a)
TL = 3000; % time interval used for the phase portraits, Figs. (b) and (c)

color_ex = [0 0 1];     % color for the original signal
color_rc = [1 0.3 0];   % color for the NGRC reconstructed signal

segLength = 50; % number of points drawn with the same color tone

% Transparency parameters for the phase portraits, Figs. (b) and (c)
a1 = 0.03;
a2 = 0.05;
 % age = k/N;
 % alpha = a1 + a2*age;   % older segments are fainter, newer segments darker

Qp = setParam(lw2, lw3, lwd, color_ex, color_rc, a1, a2, segLength, Dt, Nwind, TN, TL);

% Load data containing time, original signal, and NGRC reconstructed signal
s = load('TUV.mat');
TUV = s.TUV;

Qd = setData(TUV);

% Build the figure while hidden to reduce rendering overhead during plotting,
% especially for the segmented phase-portrait curves.
fig20 = figure(20);
set(fig20, 'Visible', 'off');

Plot_Short_Long_Intersp_map(Qp, Qd);
Plot_spectra(Qp, Qd);

set(fig20, 'Units', 'centimeters')
set(fig20, 'PaperSize', [8.5 11]);
set(fig20, 'Position', [2, 2, 8.5, 11]);

savefig(fig20, 'Fig_Mnt.fig')

saveas(gcf,'Fig_Mnt','epsc')

disp("FIN");
set(fig20, 'Visible', 'on');
% drawnow limitrate;

    return;

end % mainT_chgpt
