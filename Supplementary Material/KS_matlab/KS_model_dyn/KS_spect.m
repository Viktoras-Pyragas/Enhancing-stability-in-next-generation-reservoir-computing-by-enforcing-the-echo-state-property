%--------------------------------------------------------------------------
% Here we learn the NGRC DDE system by the Kuramoto Syvashinsky system;
% The learning can be performed with several values of Lc;
% Afterwards, we integrate the closed NGRC DDE syste for each value of
% Lc, and compute the corresponding values of Tesc (escape time)
% The learning is performed with the spectral component a1(t);
%--------------------------------------------------------------------------
%--NOTE!-------------------------------------------------------------------
% This program can be used for scanning Lc and computing corresponding
% escape times Tesc; However, we here restrict ourselves by a single value of Lc
% and compute (record) the dynamics of the NGRC;
% One can also compute (record) the NGRC dynamics without imposing
% requirement on Lc;
%--------------------------------------------------------------------------

close all
clearvars
format long;

tic;

%Ns=3;
%epsnorm=0.01; % the time series will be normalized 
epsnorm=0.0; % the time series will be normalized 
% in the range [-1+epsnorm, 1-epsnorm];
Dmx=2.0; % maximum value, criteria for escape of NGRC from solution of exact Roessler system;
%--------------------------------------------------------------------------
% With Lc:
     Lc_tot=-3.0;
    NLc=size(Lc_tot,2)
%ur=pr
%--------------------------------------------------------------------------
% Without Lc:
   % Lc_tot=0.0; % Tesc_without_Lc.mat; bet=1e-6; tpr=2000;
   % NLc=1;
%--------------------------------------------------------------------------
h=0.01; % time step; (Kesto)
tau=0.17; % from average mutual information, first minimum;
Ntau=round(tau/h); % number of steps in delay time tau
TL=1000; % Learning interval; (Kesto)
L=round(TL/h); % Number of steps in the learning interval
bet=1e-7; % ridge regression parameter; (Mano)
degree=12; % degree of nonlinearity
Ndel=4; % number of delays
m=Ndel+1; % embedding dimension
%------------------------------------------------------------
% This data file is computed and recorded in the prgram
% ks_spectral_ak_v.m
load('ks_spectral_a1.mat'); % time series og original system
%------------------------------------------------------------
X0t=a(:,1); % total array of exact signal a1(t)
Texp=0.0; % initial time shift
nTexp=round(Texp/h)
X0=X0t(1+nTexp:L+nTexp,1); % time series for learning;

Zmn=min(X0)
Zmx=max(X0)
am=2*(1-epsnorm)/(Zmx-Zmn)
Z=am*(X0.'-Zmn)-1+epsnorm; % renormalized signal x(t)
szZ=size(Z)
XX=zeros(L-Ndel*Ntau-1,Ndel+1); 
% Constructing embedding delayed space:
% [x(t-tau),x(t-2*tau),..,x(t-Ndel*tau)]
for nd=1:Ndel+1 
    XX(:,nd)=Z(1,(Ndel-(nd-1))*Ntau+1:L-(nd-1)*Ntau-1).'; % Mano
end
LL=size(XX,1)

D=(Z(1,Ndel*Ntau+2:L)-Z(1,Ndel*Ntau+1:L-1))/h; % dynamics of time derivative dx/dt
% R(1:LL,1)=ones(LL,1);
R=ones(LL,1); 
for nd=1:degree
    R=[R,generateChebyshevs(XX,nd)]; 
end
R=R.'; % This matrix will be used for optimization;
LengthW=size(R,1)
DR=zeros(LengthW,1);
% DR(1,1)=0.0 % since derivative of 1 is zero; (!!!)
sumt=1; 
%sumt=0;
for nd=1:degree
    DRnd=generateChebyshevsDiffLyap(XX,nd);
    ld=length(DRnd);
    DR(sumt+1:sumt+ld,1)=DRnd;
    sumt=sumt+ld;
end

RR=R*R.'; % Composing the matrices for optimization problem:
% With Lc requirement included:|
A=[RR, DR;
   DR.', 0];
%A=RR; % with no requirement on Lc
d=size(A,1);
A=A+bet*eye(d);
Wt=[];
for nLc=1:NLc    
    lam=Lc_tot(1,nLc);
B1=[R*D.';lam]; % with Lc requirement included|
%B1=R*D.'; % with no requirement on Lc
W1 = (A\B1).'; % the answer; row-vector;
szW1=size(W1)
disp('Lagrange multiplier')
disp(W1(end))
W=W1(1,1:end-1); % with Lc requirement included
%W=W1; % with no requirement on Lc
Wt=[Wt; W]; % with several values of Lc there may be several output vectors;
Error=sqrt(mean((D-W*R).^2))/std(D);
disp('Fitting Error')
disp(Error)
disp('Conditional Lyapunov')
disp(W*DR)
end

ND=size(D,2);
tD=h*(0:ND-1);
Xr=W*R;
figure
hold on
% Comparing the derivatives during learning; 
plot(tD,D,'b-'); % dynamics of exact derivative
plot(tD,Xr,'r-'); % dynamics of learned derivative
hold off

%% Climate replication after learning
tpr=5000; % replication interval
Tpr=(0:h:tpr);
Npr=length(Tpr); % length of replication array

%Tprad=10000;
generate_all_exp_matreces(m,degree); % generating the exponents matrices;
load('matr.mat','matr');
%Tprad_tot=(0:10:90);
% Tprad_tot=[0 10];
Tprad_tot=0;
%NTpr=size(Tprad_tot,2);
NTpr=1;
%ur=pr
Tesc_tot=zeros(NLc,NTpr);

nTprx=1;
Tesc=[];

for nLc=1:NLc
    W=Wt(nLc,:);
    nLc_v=nLc
    Lc=Lc_tot(1,nLc)
    tic;   
    countv=0;
    for nTpr=1:NTpr
        countv=countv+1;
        if countv==1
            toc;
        countv=0;
        nLc_v=nLc
        nTpr_v=nTpr
            tic;
        end

dNTpr=round(Tprad_tot(nTpr)/h);
Y1=X0t(L+1+nTexp+dNTpr:L+nTexp+dNTpr+Npr,1);
XP=zeros(1,Npr); % array of solution of NGRC DDE (will be computed below)
U=am*(Y1(:,1)-Zmn)-1+epsnorm; % exact solution (renormalized)
XP(1,1:Ndel*Ntau+1)=U(1:Ndel*Ntau+1,1).'; % initial condition for NGRC DDE system;
xx=zeros(1,m);
set1=1;

Tmx1=(Npr-Ndel*Ntau)*h;
jdx=(Ndel*Ntau+1:Npr);
Njdx=Npr-Ndel*Ntau;
jt=0;
% Integration of the NGRC DDE system:
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
    D_esc=abs(XP(1,j));
    if (D_esc>Dmx)&&(set1==1)
        Tesc=h*(j-Ndel*Ntau); % escape time
        Tesc_tot(nLc,nTpr)=Tesc;
        set1=0;
    end
    xx(1,2:Ndel+1)=XP(1,j-(1:Ndel)*Ntau); % delayed coordinates
    r=ones(LengthW,1);
    sumt=1;
    %sumt=0;
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
end % while set1==1
    if (nTpr==nTprx)&&(nLc==1)
        XPt=XP;
    end
    end % for nTpr=1:NTpr
    toc;
end % for nLc=1:NLc

Tpr1=Tpr(1,Ndel*Ntau+1:Npr).';
U1=U(Ndel*Ntau+1:Npr,1);
XP1=XP(1,Ndel*Ntau+1:Npr).';

Xtot=[Tpr1 U1 XP1];

save('Tesc_Lc.mat','Tesc_tot','Lc_tot'); % with Lc
%save('Tesc_without_Lc.mat','Tesc_tot','Lc_tot'); % without Lc

%save('XPt_record.mat','Tpr','XPt');
save('XPt_record2.mat','Xtot');

dx=abs(XP1-U1);
indx2=find(dx>0.02,1,'first');
Tvp=indx2*h;

disp('Valid prediction time:');
fprintf('%f \n',Tvp)
%Tvp=indx2*h


figure
hold on
plot(Xtot(:,1),Xtot(:,2),'r-');
plot(Xtot(:,1),Xtot(:,3),'b-');
hold off

figure
subplot(1,2,1)
plot(U(Ndel*Ntau+1:Npr,1),U((Ndel-1)*Ntau+1:Npr-Ntau,1));
ylim([-1 1])
subplot(1,2,2)
plot(XP(1,Ndel*Ntau+1:Npr),XP(1,1+(Ndel-1)*Ntau:Npr-Ntau));
ylim([-1 1])

figure
subplot(2,1,1)
hold on
plot(Tpr(1,Ndel*Ntau+1:end).',U(Ndel*Ntau+1:end,1),'-r');
plot(Tpr(1,Ndel*Ntau+1:end).',XP(1,Ndel*Ntau+1:end),'-b');
xlim([0 1000]);
ylim([-1 1])
subplot(2,1,2)
del=U(Ndel*Ntau+1:Npr,1)-XP(1,1+Ndel*Ntau:Npr).';
plot(Tpr(1,Ndel*Ntau+1:Npr).',del);


rmsed1=sqrt(mean(del.^2))

if isempty(Tesc)
    disp('No escape');
else
    disp('Escape time:');
    disp(Tesc);
end
%lamc_v=lam
 toc;
