function [PSD_exact,PSD_ngrc,f] = plot_MG_psd(Xtot)

%--------------------------------------------------
% Xtot=[Tpr1 U1 XP1]; % source array
Tpr1=Xtot(:,1);
Tpr1=Tpr1-Tpr1(1,1);
U1=Xtot(:,2);
XP1=Xtot(:,3);

% Power Spectral Density comparison using Welch algorithm
% Exact system: (Tpr1, U1)
% NGRC system:  (Tpr1, XP1)

t = Tpr1(:);
u_exact = U1(:);
u_ngrc  = XP1(:);

% Keep only samples that are valid in all three vectors
valid_idx = isfinite(t) & isfinite(u_exact) & isfinite(u_ngrc);
t = t(valid_idx);
u_exact = u_exact(valid_idx);
u_ngrc  = u_ngrc(valid_idx);

% Sampling frequency from time vector
dt = mean(diff(t));
Fs = 1 / dt;

% Remove mean before spectral analysis
u_exact = u_exact - mean(u_exact);
u_ngrc  = u_ngrc  - mean(u_ngrc);

% Welch parameters
N = length(t);
%window_length = floor(N / 8);
window_length = floor(N / 128);
window = hamming(window_length);
noverlap = floor(window_length / 2);
nfft = max(1024, 2^nextpow2(window_length));

% Compute PSDs using Welch algorithm
[PSD_exact, f] = pwelch(u_exact, window, noverlap, nfft, Fs);
[PSD_ngrc,  ~] = pwelch(u_ngrc,  window, noverlap, nfft, Fs);


end