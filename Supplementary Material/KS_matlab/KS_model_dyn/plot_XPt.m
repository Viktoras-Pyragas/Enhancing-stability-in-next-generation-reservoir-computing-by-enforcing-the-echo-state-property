h=0.01;
tpr=5000; % replication interval
%tpr=40000; % replication interval
Tpr=(0:h:tpr);
Npr=length(Tpr); % length of replication array
%--------------------------------------------------
% Tpr1=Tpr(1,Ndel*Ntau+1:Npr).';
% U1=U(Ndel*Ntau+1:Npr,1);
% XP1=XP(1,Ndel*Ntau+1:Npr).';
%Xtot=[Tpr1 U1 XP1];
%save('XPt_record2.mat','Xtot');
%load('XPt_record2_v1.mat');
% Lc_tot=-1.0; h=0.01; tau=0.12; TL=1000; bet=1e-4;
% degree=10; Ndel=4;
%--------------------------------------------------
% Tpr1=Tpr(1,Ndel*Ntau+1:Npr).';
% U1=U(Ndel*Ntau+1:Npr,1);
% XP1=XP(1,Ndel*Ntau+1:Npr).';
%Xtot=[Tpr1 U1 XP1];
%save('XPt_record2.mat','Xtot');
 load('XPt_record2.mat');
% Lc_tot=-3.0; h=0.01; tau=0.17; TL=1000; bet=1e-7;
% degree=12; Ndel=4;
% Parameters of renormalization;
%Mx=1.706986960659342;
%load('XPt_record2.mat');
epsnorm=0.0; % no shift on margins
Zmn=-1.706986960659342;
Zmx=1.704461463375478;
am=2*(1-epsnorm)/(Zmx-Zmn);
%--------------------------------------------------
Tpr1=Xtot(:,1);
U1=Xtot(:,2);
XP1=Xtot(:,3);
%-----------------------------------------------
%Error=sqrt(mean((D-W*R).^2))/std(D);
Error=sqrt(mean((U1-XP1).^2))/std(U1)

h=0.01;
mf=1;
dN=mf*50;
N0=mf*100;
N1=mf*8100;

TL=1000; % Learning interval; (Kesto)
L=round(TL/h); % Number of steps in the learning interval


Xnow=XP1(dN+1:end,1);
Xdel=XP1(1:end-dN,1);

% Xt(:,1:2)
load('ks_spectral_a1.mat');
Xt=[t(:,1) a(:,1)]; % dynamics of the first harmonic
%szXt=size(Xt)
Tt=Xt(:,1);
Xt=Xt(:,2);
Xt=am*(Xt-Zmn)-1+epsnorm; % renormalized signal x(t)

nTexp=0; % initial time shift
dNTpr=0;
X0=Xt(L+1+nTexp+dNTpr:L+nTexp+dNTpr+Npr,1);
T0=Tt(L+1+nTexp+dNTpr:L+nTexp+dNTpr+Npr,1);

Xnow_ex=Xt(dN+1:end,1);
Xdel_ex=Xt(1:end-dN,1);

Xnow=Xnow(N0:N1,1);
Xdel=Xdel(N0:N1,1);

Xnow_ex=Xnow_ex(N0:N1,1);
Xdel_ex=Xdel_ex(N0:N1,1);

figure
hold on
subplot(1,2,1);
plot(Xnow_ex,Xdel_ex,'b-');
subplot(1,2,2);
plot(Xnow,Xdel,'r-');
hold off

Tpr1=Tpr1-Tpr1(1,1);
figure
hold on
plot(Tpr1,XP1,'r--'); % time series of NGRC
plot(Tpr1,U1,'b-'); % time series of exact system
xlim([0 10]);
hold off

