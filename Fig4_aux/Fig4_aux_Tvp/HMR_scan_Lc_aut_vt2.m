%--------------------------------------------------------------------------
% Here we learn the NGRC DDE system by the HMR system;
% The learning is performed with several values of Lc;
% Afterwards, we integrate the closed NGRC DDE syste for each value of
% Lc, and compute the corresponding values of Tesc (valid prediction)
% The learning is performed with the variable x(t);
%--------------------------------------------------------------------------
close all
clearvars
format long;

tic;

Ns=3;
Dmx=0.02; % maximum value, criteria for escape of NGRC from solution of exact Roessler system;
%lam=-1; % desired conditional Lyapunov exponent; (Kesto)
%lam=-1.0; % desired conditional Lyapunov exponent; (Mano bandymui)
%--------------------------------------------------------------------------
% With Lc:
% Lc_tot=-1.0;
%  Lc_tot=-3.2;
%  Lc_tot=-3.24;
%  Lc_tot=-3.25;
%  Lc_tot=[-0.08 -1 -2 -3 -3.2 -5 -10 -20];
%   Lc_tot=[[-1 -2 -3],[-3.2 -3.24 -3.25],[-4 -5 -10]]; % Tvp_Lc_m1_m10.mat;
%  Lc_tot=[(-1:-1:-9),(-10:-10:-50)];
% Lc_tot=[(-0.5:-0.5:-5),(-6:-1:-10),(-20:-10:-50)];
% Lc_tot=(-0.5:-0.5:-5);
% Lc_tot=(-0.5:-0.5:-2.5); % Tesc_Lc_N1e3_NLc5.mat
% Lc_tot=[(-6:-1:-10),(-20:-10:-50)];
% Lc_tot=[-3.0 -3.5]; % Tesc_Lc_N1e3_NLc2.mat
% Lc_tot=[(-4.0:-0.5:-5),[-6 -7]]; % Tesc_Lc_N1e3_NLc5p.mat
% Lc_tot=[-8.0 -9.0]; % Tesc_Lc_N1e3_NLc2p.mat
% Lc_tot=(-10:-10:-50); % Tesc_Lc_N1e3_NLc5pp.mat
% Lc_tot=(-0.5:-0.5:-5);
% NLc=size(Lc_tot,2)
% ur=pr
%--------------------------------------------------------------------------
% Without Lc:
  Lc_tot=0.0;
  NLc=1;
%--------------------------------------------------------------------------
%lam=-20; % desired conditional Lyapunov exponent; (Mano)
h=0.05; % time step; (Kesto)
%h=0.015; % laiko zingsnis; (Kesto)
%h=0.01; % laiko zingsnis; (Mano)
tau=1.0; % embedding delay time
Ntau=round(tau/h); % number of steps in delay time tau
TL=2000; % Learning interval; (Kesto)
%TL=45000; % Learning interval
%TL=6000; % Learning interval
L=round(TL/h); % Number of steps in the learning interval
%bet=1e-3; % ridge regression parameter; (Kesto)
%bet=1e-7; % ridge regression parameter; (Mano)
bet=1e-6; % ridge regression parameter; (Mano)
%bet=1e-4; % ridge regression parameter; (Mano)
%bet=1e-8; % ridge regression parameter; (Mano)
%degree=9; % degree of nonlinearity
degree=8; % degree of nonlinearity
%degree=5; % degree of nonlinearity
Ndel=5; % number of delays
%Ndel=4; % number of delays
%Ndel=3; % number of delays
%Ndel=2; % number of delays
m=Ndel+1; % embedding dimension
y0=[-1.000698154309275 -4.051859895329726 3.149054785532131].'; % initial cinditions
Tin=150; % Change of initial conditions on the strange attractor
Lin=round(Tin/h);
[X0,Y0,Z0]=HMR_signal(h,Lin,y0);
y0=[X0(end),Y0(end),Z0(end)].';
[X0,Y0,Z0]=HMR_signal(h,L,y0);
%[X0,Y0,Z0,~,~,~]=Roessler_signal(h,L);
Zmn=min(X0);
Zmx=max(X0);
am=2/(Zmx-Zmn);
Z=am*(X0.'-Zmn)-1; % renormalized signal x(t)

szZ=size(Z)
XX=zeros(L-Ndel*Ntau-1,Ndel+1); % Mano
%XX=Z(1,Ndel*Ntau+1:L-1).'; % x(t) % Kesto
%szXX=size(XX)
% Constructing embedding delayed space:
% [x(t-tau),x(t-2*tau),..,x(t-Ndel*tau)]
for nd=1:Ndel+1 % Mano (1:Ndel -> Kesto)
    %XX=[XX, Z(1,(Ndel-nd)*Ntau+1:L-nd*Ntau-1).']; % Kesto
    XX(:,nd)=Z(1,(Ndel-(nd-1))*Ntau+1:L-(nd-1)*Ntau-1).'; % Mano
end
LL=size(XX,1)

D=(Z(1,Ndel*Ntau+2:L)-Z(1,Ndel*Ntau+1:L-1))/h; % dynamics of time derivative dx/dt

yend=[X0(end) Y0(end) Z0(end)]; % state of system at the end of learning

R(1:LL,1)=ones(LL,1);
R=ones(LL,1); % Kesto
for nd=1:degree
    R=[R,generateChebyshevs(XX,nd)]; % Kesto 
    %R=[R,generateChebyshevs_v(XX,nd)]; % Kesto    
end
R=R.'; % This matrix will be used for optimization;
LengthW=size(R,1)

DR=zeros(LengthW,1);
% DR(1,1)=0.0 % since derivative of 1 is zero; (!!!)
sumt=1; % Kesto
for nd=1:degree
    DRnd=generateChebyshevsDiffLyap(XX,nd);
    %DRnd=generateChebyshevsDiffLyap_v(XX,nd);
    ld=length(DRnd);
    DR(sumt+1:sumt+ld,1)=DRnd;
    sumt=sumt+ld;
end

RR=R*R.'; % Composing the matrices for optimization problem:
% A=[RR,DR];
% A=[A;[DR.',0]];
% A=[RR, DR;
%    DR.', 0];
A=RR;
d=size(A,1);
A=A+bet*eye(d);
Wt=[];
for nLc=1:NLc
    lam=Lc_tot(1,nLc);
%B1=[R*D.';lam];
B1=R*D.';
W1 = (A\B1).'; % the answer; row-vector;
szW1=size(W1)
%W1 = (pinv(A)*B1).';
disp('Lagrange multiplier')
disp(W1(end))
%W=W1(1,1:end-1);
W=W1;
Wt=[Wt; W];
Error=sqrt(mean((D-W*R).^2))/std(D);
disp('Fitting Error')
disp(Error)
disp('Conditional Lyapunov')
disp(W*DR)
end
%disp(W)
%T=linspace(0,hT*(NT-1),NT);
% ur=pr
%% Climate replication after learning
tpr=1000; % replication interval
%tpr=1000; % replication interval
%tpr=10000; % replication interval
Tpr=(0:h:tpr);
Npr=length(Tpr); % length of replication array

r=0.006;
nu=4;
kap=-1.6;
I=3.2;

P.r=r;
P.nu=nu;
P.kap=kap;
P.I=I;

% a=0.2; 
% b=0.2;
% %b=0.6;
% c=5.7;
% P.a=a; P.b=b; P.c=c;
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);

%Tprad=10000;
generate_all_exp_matreces(m,degree); % generating the exponents matrices;
load('matr.mat','matr');
%Tprad_tot=(1000:100:1900); % for HMR
%Tprad_tot=(1000:500:5500); % for HMR
%Tprad_tot=(1000:500:10500); % for HMR
%Tprad_tot=(1000:500:50500); % for HMR
%Tprad_tot=(1000:500:500500); % for HMR
%Tprad_tot=(10000:1000:19000);
%Tprad_tot=(10000:1000:29000);
%Tprad_tot=(10000:1000:109000);
%Tprad_tot=(10000:1000:109000);
%Tprad_tot=(10000:1000:1009000);
% Tprad_tot=(10000:1000:109000); % n1
% Tprad_tot=(110000:1000:209000); % n2
% Tprad_tot=(210000:1000:309000); % n3
% Tprad_tot=(310000:1000:409000); % n4
% Tprad_tot=(410000:1000:509000); % n5
  Tprad_tot=(10000:1000:509000); % NTpr=500;
NTpr=size(Tprad_tot,2);
%ur=pr
Tesc_tot=zeros(NLc,NTpr);

y0pr=zeros(NTpr,Ns);
for nTpr=1:NTpr
if nTpr==1        
y0=yend;
Tprad=Tprad_tot(1,nTpr)
else
Tprad=Tprad_tot(1,nTpr)-Tprad_tot(1,nTpr-1);    
end
[~,Y1] = ode45(@(t,y)sistema(t,y,P),[0 Tprad],y0,opts); % warming up
y0=Y1(end,:);
y0pr(nTpr,:)=y0;
end

for nLc=1:NLc
    W=Wt(nLc,:);
    nLc_v=nLc
    Lc=Lc_tot(1,nLc)
    tic;   
    countv=0;
    for nTpr=1:NTpr
        countv=countv+1;
        if countv==10
            toc;
        countv=0;
        nLc_v=nLc
        nTpr_v=nTpr
            tic;
        end
y0=y0pr(nTpr,:);
[~,Y1] = ode45(@(t,y) sistema(t,y,P), Tpr,y0,opts);

XP=zeros(1,Npr); % array of solution of NGRC DDE (will be computed below)
U=am*(Y1(:,1)-Zmn)-1; % exact solution (renormalized)

%XP(1,(Ndel-1)*Ntau+1)=U(Ndel*Ntau+1,1); % initial condition for NGRC DDE system;
XP(1,1:Ndel*Ntau+1)=U(1:Ndel*Ntau+1,1).'; % initial condition for NGRC DDE system;
xx=zeros(1,m);
set1=1;

jdx=(Ndel*Ntau+1:Npr-1);
Njdx=Npr-1-Ndel*Ntau;
jt=0;
% Integration of the NGRC DDE system:
%for j=(Ndel-1)*Ntau+1:Npr-1
%for j=Ndel*Ntau+1:Npr-1
while set1==1
    jt=jt+1;
    if (set1==1)
        if jt<Njdx
            j=jdx(1,jt);
        else
            j=jdx(1,end);
            set1=0;
        end
    end    
    xx(1,1)=XP(1,j); % current state
    D_esc=abs(XP(1,j)-U(j,1));
    if (D_esc>Dmx)&&(set1==1)
        Tesc=h*(j-Ndel*Ntau); % escape time
        Tesc_tot(nLc,nTpr)=Tesc;
        set1=0;
    end
    %xx(1,2:Ndel+1)=U(j-(0:Ndel-1)*Ntau,1).';
    xx(1,2:Ndel+1)=XP(1,j-(1:Ndel)*Ntau); % delayed coordinates
    r=ones(LengthW,1);
    sumt=1;
    for nd=1:degree
    matrnd=matr(nd).f;
    rn=generateCheb_feat_vect(xx,matrnd,nd);
    ln=length(rn);
    r(sumt+1:sumt+ln,1)=rn;
    sumt=sumt+ln;
    end
    if set1==1
    XP(1,j+1)=XP(1,j)+h*(W*r);
    end
    % Restriction:
    % if(abs(XP(1,j+1))>1)
    % XP(1,j+1)=sign(XP(1,j));       
    % end
end % while set1==1
    
    end % for nTpr=1:NTpr
    toc;
end % for nLc=1:NLc

Tvp_tot=Tesc_tot;
%save('Tvp_Lc.mat','Tvp_tot','Lc_tot'); % with Lc
save('Tvp_without_Lc.mat','Tvp_tot','Lc_tot'); % without Lc

figure
subplot(1,2,1)
%plot(U((Ndel-1)*Ntau+1:Npr,1),XP(1,(Ndel-1)*Ntau+1:Npr).')
plot(U(Ndel*Ntau+1:Npr,1),U((Ndel-1)*Ntau+1:Npr-Ntau,1));
ylim([-1 1])
subplot(1,2,2)
%plot(U(Ndel*Ntau+1:Npr,1),XP(1,1+(Ndel-1)*Ntau:Npr-Ntau).')
plot(XP(1,Ndel*Ntau+1:Npr),XP(1,1+(Ndel-1)*Ntau:Npr-Ntau));
ylim([-1 1])

figure
subplot(2,1,1)
hold on
% plot(Tpr(1,(Ndel-1)*Ntau+1:end).',U((Ndel-1)*Ntau+1:end,1),'-b')
% plot(Tpr(1,(Ndel-1)*Ntau+1:end).',XP(1,(Ndel-1)*Ntau+1:end),'-r')
plot(Tpr(1,Ndel*Ntau+1:end).',U(Ndel*Ntau+1:end,1),'-b');
plot(Tpr(1,Ndel*Ntau+1:end).',XP(1,Ndel*Ntau+1:end),'-r');
%xlim([0 200]);
xlim([0 1000]);
ylim([-1 1])
subplot(2,1,2)
%del=U(Ndel*Ntau+1:Npr,1)-XP(1,1+(Ndel-1)*Ntau:Npr-Ntau).';
del=U(Ndel*Ntau+1:Npr,1)-XP(1,1+Ndel*Ntau:Npr).';
plot(Tpr(1,Ndel*Ntau+1:Npr).',del);


rmsed1=sqrt(mean(del.^2))
disp('Escape time:');
disp(Tesc);
lamc_v=lam
 toc;
 
function yp = sistema(~,y,P)
% Hindmarsh Rose system
  yp = [y(2)-y(1)^3+3*y(1)^2-y(3)+P.I;
        1-5*y(1)^2-y(2);
        P.r*(P.nu*(y(1)-P.kap)-y(3))];
end
