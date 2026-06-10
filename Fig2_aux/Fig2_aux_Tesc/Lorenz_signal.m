function  [X,Y,Z]=Lorenz_signal(hT,NT,y0)
% Generating time series of X,Y,Z of the Roesslerio system
% hT time step
% NT number of points (length of the signal)
% y0 intial condition 
P.sig=10; 
P.r=28;
P.b=8/3; % Lorenz parameters
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
T=linspace(0,hT*(NT-1),NT);
[~,Y1] = ode45(@(t,y) sistema(t,y,P), T,y0,opts);
%format long
disp(Y1(end,:))
%format short
X=Y1(:,1);
Y=Y1(:,2);
Z=Y1(:,3);
end

function yp = sistema(~,y,P)
% Lorenz sistema
yp=[-P.sig*(y(1)-y(2));
     P.r*y(1)-y(2)-y(1)*y(3);
     y(1)*y(2)-P.b*y(3)];
end