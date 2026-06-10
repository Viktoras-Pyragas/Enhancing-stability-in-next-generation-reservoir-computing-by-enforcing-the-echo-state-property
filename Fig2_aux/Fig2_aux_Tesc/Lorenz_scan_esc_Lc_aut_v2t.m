%--------------------------------------------------------------------------
% Here we learn the NGRC DDE system by the Lorenz system;
% The learning is performed with several values of Lc;
% Afterwards, we integrate the closed NGRC DDE syste for each value of
% Lc, and compute the corresponding values of Tesc (escape time from attractor)
% The learning is performed with the variable x(t);
%--------------------------------------------------------------------------
close all
clearvars
format long;

tic;

Ns=3;
%Dmx=0.02; % maximum value, criteria for escape of NGRC from solution of exact Roessler system;
Dmx_esc=2.0;
%lam=-1; % desired conditional Lyapunov exponent; (Kesto)
%lam=-1.0; % desired conditional Lyapunov exponent; (Mano bandymui)
%--------------------------------------------------------------------------
% With Lc:
%Lc_tot=-10.0;
Lc_tot=-9.0;
%Lc_tot=-5.0;
%Lc_tot=[(-7:-1:-9)];
%Lc_tot=[(-0.1:-0.1:-0.4),(-0.5:-0.5:-2.5),(-3:-1:-5),(-10:-10:-50),-100];
%Lc_tot=[(-0.1:-0.1:-0.4)];
%Lc_tot=[(-0.5:-0.5:-2.5)];
%Lc_tot=[(-3.0:-0.5:-5)];
%Lc_tot=[(-4.6:-0.1:-4.9)];
%Lc_tot=[(-0.5:-0.5:-5)];
%Lc_tot=[(-0.5:-0.5:-5),(-6:-1:-10),(-20:-10:-50),(-100:-50:-200)];
%Lc_tot=(-0.5:-0.5:-5);
%Lc_tot=(-0.5:-0.5:-3);
%Lc_tot=[(-3.5:-0.5:-5),(-6:-1:-10),(-20:-10:-50),(-100:-50:-200)];
% Lc_tot=[(-6:-1:-10),(-20:-10:-50)];
%Lc_tot=(-100:-50:-200);
NLc=size(Lc_tot,2)
%ur=pr
%--------------------------------------------------------------------------
% Without Lc:
 % Lc_tot=0.0;
 % NLc=1;
%--------------------------------------------------------------------------
%lam=-20; % desired conditional Lyapunov exponent; (Mano)
h=0.01; % time step; (I&K)
%h=0.015; % laiko zingsnis; (Kesto)
%h=0.01; % laiko zingsnis; (Mano)
tau=0.15; % embedding delay time
Ntau=round(tau/h); % number of steps in delay time tau
TL=300; % Learning interval; (Kesto)
%TL=45000; % Learning interval
%TL=6000; % Learning interval
L=round(TL/h); % Number of steps in the learning interval
bet=1e-7; % ridge regression parameter; (Mano)
%bet=1e-8; % ridge regression parameter; (Mano)
%bet=1e-6; % ridge regression parameter; (Mano)
degree=7; % degree of nonlinearity
%degree=9; % degree of nonlinearity
Ndel=5; % number of delays
%Ndel=7; % number of delays
m=Ndel+1; % embedding dimension
y0=[-6.924443285060142 -12.968332920785745  10.820682533915615].'; % initial conditions
[X0,Y0,Z0]=Lorenz_signal(h,L,y0);
M=19.3667  % max computed from a very long (T=10000) time interval;
Z=X0.';
Z=Z/M;

szZ=size(Z)
XX=zeros(L-Ndel*Ntau-1,Ndel+1); % Mano
%XX=Z(1,Ndel*Ntau+1:L-1).'; % y(t) % Kesto
%szXX=size(XX)
% Constructing embedding delayed space:
% [y(t-tau),y(t-2*tau),..,y(t-Ndel*tau)]
for nd=1:Ndel+1 % Mano (1:Ndel -> Kesto)
    %XX=[XX, Z(1,(Ndel-nd)*Ntau+1:L-nd*Ntau-1).']; % Kesto
    XX(:,nd)=Z(1,(Ndel-(nd-1))*Ntau+1:L-(nd-1)*Ntau-1).'; % Mano
end
LL=size(XX,1)

D=(Z(1,Ndel*Ntau+2:L)-Z(1,Ndel*Ntau+1:L-1))/h; % dynamics of time derivative dy/dt

yend=[X0(end) Y0(end) Z0(end)]; % state of system at the end of learning

R(1:LL,1)=ones(LL,1);
R=ones(LL,1); % Kesto
for nd=1:2:degree
    R=[R,generateChebyshevs(XX,nd)]; % Kesto    
end
R=R.'; % This matrix will be used for optimization;
LengthW=size(R,1)

DR=zeros(LengthW,1);
% DR(1,1)=0.0 % since derivative of 1 is zero; (!!!)
sumt=1; % Kesto
for nd=1:2:degree
    DRnd=generateChebyshevsDiffLyap(XX,nd);
    ld=length(DRnd);
    DR(sumt+1:sumt+ld,1)=DRnd;
    sumt=sumt+ld;
end

RR=R*R.'; % Composing the matrices for optimization problem:
% A=[RR,DR];
% A=[A;[DR.',0]];
A=[RR, DR;
   DR.', 0];
%A=RR;
d=size(A,1);
A=A+bet*eye(d);
Wt=[];
for nLc=1:NLc
    lam=Lc_tot(1,nLc);
B1=[R*D.'; lam];
%B1=R*D.';
W1 = (A\B1).'; % the answer; row-vector;
szW1=size(W1)
%W1 = (pinv(A)*B1).';
disp('Lagrange multiplier')
disp(W1(end))
W=W1(1,1:end-1);
%W=W1;
Wt=[Wt; W];
Error=sqrt(mean((D-W*R).^2))/std(D);
disp('Fitting Error')
disp(Error)
disp('Conditional Lyapunov')
disp(W*DR)
end
%disp(W)
%T=linspace(0,hT*(NT-1),NT);

%% Climate replication after learning
tpr=1000; % replication interval
%tpr=10000; % replication interval
%tpr=50000; % replication interval
%tpr=500000; % replication interval
%tpr=100000; % replication interval
%tpr=25000; % replication interval
Tpr=(0:h:tpr);
Npr=length(Tpr); % length of replication array
Tesc_mx=h*(Npr-Ndel*Ntau); % maximal possible escape time;

sig=10; 
r=28;
b=8/3; % Lorenz parameters
P.sig=sig; P.r=r; P.b=b;
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);

%Tprad=10000;
generate_all_exp_matreces(m,degree); % generating the exponents matrices;
load("matr.mat","matr");
Tprad_tot=1000;
%Tprad_tot=(1000:100:1100);
%Tprad_tot=(1000:100:1900);
%Tprad_tot=(1000:100:2900);
%Tprad_tot=(1000:100:10900);
%Tprad_tot=(1000:100:100900);
%Tprad_tot=(10000:1000:109000);
%Tprad_tot=(10000:1000:1009000);
NTpr=size(Tprad_tot,2);
%Tesc_tot=zeros(NLc,NTpr);
Tesc_tot=Tesc_mx*ones(NLc,NTpr);

y0pr=zeros(NTpr,Ns);
for nTpr=1:NTpr
if nTpr==1        
y0=yend;
Tprad=Tprad_tot(1,nTpr)
else
Tprad=Tprad_tot(1,nTpr)-Tprad_tot(1,nTpr-1);    
end
[~,Y1] = ode45(@(t,y) sistema(t,y,P), [0 Tprad],y0,opts); % warming up
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
Tpr1=h*(0:Ndel*Ntau+1);
Npr1=size(Tpr1,2);
[~,Y1] = ode45(@(t,y) sistema(t,y,P), Tpr1,y0,opts);

XP=zeros(1,Npr); % array of solution of NGRC DDE (will be computed below)
U=Y1(:,1)/M; % exact solution (renormalized)

%XP(1,(Ndel-1)*Ntau+1)=U(Ndel*Ntau+1,1); % initial condition for NGRC DDE system;
XP(1,1:Ndel*Ntau+1)=U(1:Ndel*Ntau+1,1).'; % initial condition for NGRC DDE system;
xx=zeros(1,m);
set1=1;
jdx=(Ndel*Ntau+1:Npr-1);
szjdx=size(jdx,2);
% Integration of the NGRC DDE system:
%for j=(Ndel-1)*Ntau+1:Npr-1
jt=0;
%for j=Ndel*Ntau+1:Npr-1
while set1==1
    jt=jt+1;
    if (jt>szjdx)&&(set1==1)
        set1=0;
        j=jdx(1,end);
    else
        j=jdx(1,jt);
    end
    xx(1,1)=XP(1,j); % current state
    D_esc=abs(XP(1,j));
    if (D_esc>Dmx_esc)&&(set1==1)
        Tesc=h*(j-Ndel*Ntau); % escape time
        Tesc_tot(nLc,nTpr)=Tesc;
        set1=0;
    end
    %xx(1,2:Ndel+1)=U(j-(0:Ndel-1)*Ntau,1).';
    xx(1,2:Ndel+1)=XP(1,j-(1:Ndel)*Ntau); % delayed coordinates
    r=ones(LengthW,1);
    sumt=1;
    for nd=1:2:degree
    matrnd=matr(nd).f;
    rn=generateCheb_feat_vect(xx,matrnd,nd);
    ln=length(rn);
    r(sumt+1:sumt+ln,1)=rn;
    sumt=sumt+ln;
    end
    % szW=size(W)
    % szr=size(r)
    XP(1,j+1)=XP(1,j)+h*(W*r);
    % Restriction:
    % if(abs(XP(1,j+1))>1)
    % XP(1,j+1)=sign(XP(1,j));       
    % end
end
    if (nLc==1)&&(nTpr==1)
        NXP=size(XP,2);
        Tt=h*(0:NXP-1).';
        Xt=[Tt XP.'];
    end

    end % for nTpr=1:NTpr
    toc;
    
end % for nLc=1:NLc


%save('Tesc_XP.mat','Xt'); % tpr=20000;
%save('Tesc_XP2.mat','Xt'); % tpr=50000;
%save('Tesc_XP3.mat','Xt'); % tpr=100000;
%save('Tesc_XP4.mat','Xt'); % tpr=5e5;
save('Tesc_XP5.mat','Xt'); % tpr=5e5; Lc=-9.0;
%save('Tesc_attr_Lc.mat','Tesc_tot','Lc_tot'); % with Lc
%save('Tesc_attr_without_Lc.mat','Tesc_tot','Lc_tot'); % without Lc
figure
subplot(1,2,1)
%plot(U((Ndel-1)*Ntau+1:Npr,1),XP(1,(Ndel-1)*Ntau+1:Npr).')
plot(U(Ndel*Ntau+1:Npr1,1),U((Ndel-1)*Ntau+1:Npr1-Ntau,1));
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
%plot(Tpr(1,Ndel*Ntau+1:Npr1).',U(Ndel*Ntau+1:Npr1,1),'-b');
plot(Tpr(1,Ndel*Ntau+1:end).',XP(1,Ndel*Ntau+1:end),'-r');
%xlim([0 200]);
xlim([0 1000]);
ylim([-1 1])
subplot(2,1,2)
%del=U(Ndel*Ntau+1:Npr,1)-XP(1,1+(Ndel-1)*Ntau:Npr-Ntau).';
del=U(Ndel*Ntau+1:Npr1,1)-XP(1,1+Ndel*Ntau:Npr1).';
plot(Tpr(1,Ndel*Ntau+1:Npr1).',del);


rmsed1=sqrt(mean(del.^2))
disp('Escape time:');
disp(Tesc);
lamc_v=lam
 toc;
 
function yp = sistema(~,y,P)
% Lorenz system
sig=P.sig; 
r=P.r;
b=P.b;
yp=[-sig*(y(1)-y(2));
     r*y(1)-y(2)-y(1)*y(3);
     y(1)*y(2)-b*y(3)];
end
