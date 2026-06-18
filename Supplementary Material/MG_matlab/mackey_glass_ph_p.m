% Chaotic Mackey-Glass system phase portrait: (P(t), P(t - tau))
clear; close all; clc;

% Model parameters for the Mackey-Glass delay differential equation.
gamma = 1.0;
beta0 = 2.0;
n = 10.0;
tau = 2.0;

% Simulation settings. A long time span lets the transient settle.
t0 = 0;
%tf = 2000;
tf = 20000;
dt = 0.01;
t_eval = t0:dt:tf;

% Constant positive history for t <= 0.
history = @(t) 1.2;

% Right-hand side of the DDE. Z(1) represents the delayed value P(t - tau).
mg_rhs = @(t, P, Z) beta0 * Z(1) / (1 + Z(1)^n) ...
                   - gamma * P(1);

tic;
disp('Integrating DDE:');
% Solve the delay differential equation over the selected time interval.
sol = dde23(mg_rhs, tau, history, [t0 tf]);

% Evaluate the numerical solution on the uniform time grid.
P = deval(sol, t_eval);

disp('Computing delay:');

% Remove the initial transient and align P(t) with P(t - tau).
transient_time = 500;
Ntau=round(tau/dt);
Ntrans=round(transient_time/dt);
P=P(1,Ntrans+1:end);
P_delay=P(1,1:end-Ntau);
P=P(1,1+Ntau:end);
t_eval=t_eval(1,1+Ntrans:end-Ntau);

toc;

% Plot the phase portrait using the current and delayed signal values.
figure('Color', 'w');
plot(P, P_delay, 'b-', 'LineWidth', 1.0);
grid on;
box on;
xlabel('P(t)');
ylabel('P(t - \tau)');
title('Mackey-Glass Phase Portrait');

% Plot the time evolution after removing the transient part.
figure('Color', 'w');
plot(t_eval, P, 'r-', 'LineWidth', 1.0);
grid on;
box on;
xlabel('t');
ylabel('P(t)');
title('Mackey-Glass Dynamics');

% Store the simulated time series as two columns: time and P(t).
Pt=[t_eval.', P.'];
save('MG_dyn.mat','Pt');
%load('MG_dyn.mat');

% Reload the saved matrix variables from memory for an additional check plot.
t_eval1=Pt(:,1);
P1=Pt(:,2);

% Display the saved trajectory data in a separate figure.
figure
hold on
plot(t_eval1,P1,'b-');
hold off
