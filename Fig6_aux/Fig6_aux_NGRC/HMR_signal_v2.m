function  [X,Y,Z]=HMR_signal_v2(hT,NT,y0,parv)
% Generuojami X,Y,Z Roesslerio signalai ilgio NT su zingsniu hT
% Signalai pernormuojami taip kad ju vidurkiai 0 o variancai lygus 1
% hT=0.05; % laiko zingsnis
% NT=10000; % Tasku skaicius (masivo ilgis) 
r=0.006;
nu=4;
kap=-1.6;

P.nu=nu;
P.r=r; 
P.kap=kap;
P.I=parv;
%Tprad=500; % Warming up time
%Tprad=2000; % Warming up time
Tprad=1000; % Warming up time
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
%y0=[-6.877777769694082   0.195379905513923   0.015978401887531]; % initial condition
%[~,Y1] = ode45(@sistema,[0 Tprad],y0,opts);
[~,Y1] = ode45(@(t,y)sistema(t,y,P), [0 Tprad],y0,opts); % warming up
y00=Y1(end,:);

[~,Y1] = ode45(@(t,y)sistema(t,y,P), [0:hT:Tprad],y00,opts); % warming up
y00=Y1(end,:);
Zmn=min(Y1(:,2));
Zmx=max(Y1(:,2));
T=linspace(0,hT*(NT-1),NT); % interval of integration
%[~,Y1] = ode45(@sistema,T,y00,opts);
[~,Y1] = ode45(@(t,y)sistema(t,y,P), T,y00,opts); % integration
X=(Y1(:,1)).';
Y=(Y1(:,2)).';
Z=(Y1(:,3)).';
% Zmn=min(Y1(:,2));
% Zmx=max(Y1(:,2));
% DX=-X-Z;
% DY=X+a*Y;
% DZ=b+Z.*(X-c);

end
%---------------------------------------------------------------------
      
function dy=sistema(~,y,P)
% HMR system;
dy=[y(2)-y(1)^3+3*y(1)^2-y(3)+P.I;
    1-5*y(1)^2-y(2);
    P.r*(P.nu*(y(1)-P.kap)-y(3))];
end

   