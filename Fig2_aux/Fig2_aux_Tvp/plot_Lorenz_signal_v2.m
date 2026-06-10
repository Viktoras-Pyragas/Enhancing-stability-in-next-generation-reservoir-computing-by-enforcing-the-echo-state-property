% Plot Lorenz signal
close all
clearvars

format long;

hT=0.01;
%hT=0.005;
%NT=100000;
%NT=200000;
%NT=10000000;
NT=50000000;
y0=[-1.585199705556874  -2.796064168194566  12.088358657189346]'; % initial cinditions
[X,Y,Z]=Lorenz_signal(hT,NT,y0);
T=linspace(0,hT*(NT-1),NT);
figure
plot(T,X)
xlabel('time') 
ylabel('x') 
title('Lorenz')


Zt=[T.' Z];
save('Uvt.mat','Zt');

Xt=[T.' X];
%save('Xvt.mat','Xt'); % h=0.01;
%save('Xvt2.mat','Xt'); % h=0.005;

%save('Xvt3.mat','Xt'); % h=0.01; NT=1e7;
%save('Xvt4.mat','Xt'); % h=0.01; NT=5e7;
save('Xvt5.mat','Xt'); % h=0.01; NT=5e7;