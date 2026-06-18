% Integrate the spectral Kuramoto-Sivashinsky equations for a_k, k = 1,...,N,
% and plot the dynamics of the first mode a_1.

clear;
clc;
close all;

N = 32;
nu = 0.0215;
dt=0.01;

tspan0=(0:dt:300);
tspan = (0:dt:10000);

% Initial condition for the spectral coefficients a_1,...,a_N.
% Change this vector if a different trajectory is needed.
a0 = zeros(N, 1);
a0(1) = 0.8;
a0(2) = 0.4;
a0(3) = 0.2;

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

fprintf('Integrating spectral KS system with N = %d and nu = %.4f...\n', N, nu);
tic;
[t, a] = ode15s(@(t, a) ksSpectralRhs(t, a, N, nu), tspan0, a0, options);
a0=a(end,:);
[t, a] = ode15s(@(t, a) ksSpectralRhs(t, a, N, nu), tspan, a0, options);
elapsedTime = toc;
fprintf('Elapsed wall time: %.3f seconds\n', elapsedTime);

figure;
plot(t, a(:, 1), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Time');
ylabel('a_1');
title('Dynamics of the First Spectral Mode a_1');

save('ks_spectral_a1.mat', 't', 'a', 'N', 'nu');

function da = ksSpectralRhs(~, a, N, nu)
    da = zeros(N, 1);

    for k = 1:N
        sum1 = 0.0;
        sum2 = 0.0;
        sum3 = 0.0;

        % sum from m = k - N to -1 of a_{-m} a_{k-m}
        if k - N <= -1
            for m = (k - N):-1
                sum1 = sum1 + a(-m) * a(k - m);
            end
        end

        % sum from m = 1 to k - 1 of a_m a_{k-m}
        if k >= 2
            for m = 1:(k - 1)
                sum2 = sum2 + a(m) * a(k - m);
            end
        end

        % sum from m = k + 1 to N of a_m a_{m-k}
        if k <= N - 1
            for m = (k + 1):N
                sum3 = sum3 + a(m) * a(m - k);
            end
        end

        da(k) = (k^2 - nu * k^4) * a(k) + 0.5 * k * (sum1 - sum2 + sum3);
    end
end
