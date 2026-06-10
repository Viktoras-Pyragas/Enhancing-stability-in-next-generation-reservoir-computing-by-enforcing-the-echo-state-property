% Plot Lorenz signal
close all
clearvars
hT=0.01;
NT=10000;
y0=[-1.585199705556874  -2.796064168194566  12.088358657189346]'; % initial cinditions
[X,Y,Z]=Lorenz_signal(hT,NT,y0);
T=linspace(0,hT*(NT-1),NT);
figure
plot(T,X)
xlabel('time') 
ylabel('x') 
title('Lorenz')
