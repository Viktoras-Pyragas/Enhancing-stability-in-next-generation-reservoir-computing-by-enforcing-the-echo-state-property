% Average Mutual Information and first minimum
% Input: scalar time series x

% Reset the MATLAB workspace, command window, and any open figures so the
% script starts from a clean state.
clear; clc; close all

% KS model
% Load the Kuramoto-Sivashinsky spectral dataset. The file provides the time
% vector t and the matrix a, whose first column is used as the scalar signal.
load('ks_spectral_a1.mat');
T=t;
x=a(:,1);

% Plot the raw scalar time series before normalization.
figure
plot(T,x)
xlabel('Time')
ylabel('X')

% Force x to be a column vector, then standardize it to zero mean and unit
% standard deviation. This makes the AMI calculation scale-independent.
x = x(:);
x = x - mean(x);
x = x / std(x);

% maxTau controls the largest delay tested. nBins controls the resolution of
% the histogram-based probability estimates used for mutual information.
maxTau = 200;          % maximum delay to test
nBins = 32;            % histogram bins

% Preallocate storage for AMI values, one value for each tested delay.
AMI = zeros(maxTau,1);

% Compare the signal with delayed copies of itself. For each delay tau, x1
% contains the original samples and x2 contains the same signal shifted by tau.
for tau = 1:maxTau
    x1 = x(1:end-tau);
    x2 = x(1+tau:end);

    % Estimate the mutual information between x(t) and x(t+tau).
    AMI(tau) = mutual_information_hist(x1,x2,nBins);
end

% Find first local minimum
% The first local minimum is commonly used as the embedding delay for
% delay-coordinate reconstruction.
tau_embed = NaN;
for tau = 2:maxTau-1
    if AMI(tau) < AMI(tau-1) && AMI(tau) < AMI(tau+1)
        tau_embed = tau;
        break
    end
end

% Plot the AMI curve. If a first local minimum was found, mark it with a red
% circle and vertical dashed line.
figure
plot(1:maxTau,AMI,'LineWidth',1.5)
hold on
if ~isnan(tau_embed)
    plot(tau_embed,AMI(tau_embed),'ro','MarkerSize',8,'LineWidth',2)
    xline(tau_embed,'r--')
    title(sprintf('Average Mutual Information, first minimum = %d',tau_embed))
else
    title('Average Mutual Information: no local minimum found')
end
xlabel('Delay \tau')
ylabel('AMI(\tau)')
grid on

% Print the selected embedding delay in the command window.
fprintf('Embedding time delay from first AMI minimum: tau = %d\n',tau_embed)

% Estimate mutual information using histogram-based probability estimates.
% x and y are paired samples, and nBins sets the number of bins per axis.
function I = mutual_information_hist(x,y,nBins)
    % Build equally spaced histogram bin edges for both variables.
    edgesX = linspace(min(x),max(x),nBins+1);
    edgesY = linspace(min(y),max(y),nBins+1);

    % Estimate marginal probabilities p(x) and p(y).
    px = histcounts(x,edgesX,'Normalization','probability');
    py = histcounts(y,edgesY,'Normalization','probability');

    % Estimate the joint probability p(x,y).
    pxy = histcounts2(x,y,edgesX,edgesY,'Normalization','probability');

    % Accumulate I = sum p(x,y) * log(p(x,y)/(p(x)*p(y))).
    I = 0;

    for i = 1:nBins
        for j = 1:nBins
            % Skip empty bins to avoid log(0) and division by zero.
            if pxy(i,j) > 0 && px(i) > 0 && py(j) > 0
                I = I + pxy(i,j)*log(pxy(i,j)/(px(i)*py(j)));
            end
        end
    end
end
%tau_embed` is the embedding delay. The AMI is computed in natural units; use `log2` instead of `log` for bits.
