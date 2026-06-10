
%------------------------------------------------------------------------
% Here we learn the NGRC DDE by the x(t) signal from the HMR neuron sytem,
% and draw the bifurcation diagram dependent on the parameter I;
% The NGRC DDE is learned on nP values of the parameter I;
% The requirement for conditional exponent Lc is also imposed for
% each value of the parameter I;
%------------------------------------------------------------------------
close all
clearvars
format long;

Ns=3; % order of the original system
%nP=4; % number of parameters for learning;
nP=5; % number of parameters for learning;
Nb=101; % number of par-ters for computing the bifurcation
%lam=-ones(1,nP)*2; %[-10 -10]; % size(lam,2)=nP; % Mano
lam=-ones(1,nP)*1; %[-10 -10]; % size(lam,2)=nP; % Kesto;
h=0.05; % time step; (Kesto)
tau=1.0; % embedding delay time
Ntau=round(tau/h); % number of steps in delay time tau
TL=2000; % Learning interval % Kesto
%TL=300; % Learning interval % Mano
L=round(TL/h); % Number of steps in the learning interval
%bet=1e-3; % ridge regression parameter % Kesto;
%bet=1e-4;
bet=1e-6; % ridge regression parameter % (Kesto);
%bet=1e-7;
%degree=6;
%degree=4; % degree of nonlinearity; % I&K;
degree=8; % degree of nonlinearity; % Kesto
%Ndel=3; % number of delays
%Ndel=6; % number
Ndel=5; % number of delays
m=Ndel+1; % embedding dimension
%pbv=linspace(0.6,0.65,nP); % set of parameter values for learning; % I&K
%pbv=linspace(0.6,0.95,nP);
%pbv=linspace(0.2,1.0,nP);
Ivt=linspace(3.0,3.5,nP); % set of parameter values for learning; % (Mano)
%pbv=[0.6 0.8 0.9 0.95];
%pM=linspace(0.2,1.0,Nb); % array of par-ter b to compute the bif-tion;
pM=linspace(3.0,3.5,Nb); % array of par-ter b to compute the bif-tion;
Imx=pM(1,end); % maximal value; will be used for renormalization;
%-----------------------------------------------------
% bmd=[3.0 3.5]; % edges of the parameter (I) variation; 
% bm=2/(bmd(2)-bmd(1)); % normalization constant for parameter I;
%Iz=Ivt;
y0=[-1.000698154309275  -4.051859895329726   3.149054785532131].'; % initial cinditions
Tprad=1000;
disp('Initial computation of max and min:');
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
r=0.006;
nu=4;
kap=-1.6;

P.I=3.0;
P.nu=nu;
P.r=r; 
P.kap=kap;

y0P=zeros(nP,Ns); % initial conditions for original system;
y1=y0;
    tic;
    %P.b=bmd(1);
    [~,Y1] = ode45(@(t,y)sistema(t,y,P),[0 Tprad],y1,opts); % warming up
    y1=Y1(end,:);
    Tpwarm=(0:h:Tprad);
    [~,Y1] = ode45(@(t,y)sistema(t,y,P),Tpwarm,y1,opts); % warming up
    y1=Y1(end,:);
    % Finding the ranges for z(t):
    Zmn=min(Y1(:,1))
    Zmx=max(Y1(:,1)) % exact maximum;
%-------------------------------------------------------------------------    
%     Zmx=48.0; % changed manualy;
%-------------------------------------------------------------------------    
    am=2/(Zmx-Zmn)
 %   ur=pr
    toc;
    y0=y1; % the initial conditions with I=3.0;
     %y1=y0;
disp('Computing initial conditions for learning:');
for np=1:nP
    tic;
      np_v0=np
     Iv=Ivt(1,np); % setting parameter value for learning;
     P.I=Iv;
     %y1=y0;
    [~,Y1] = ode45(@(t,y)sistema(t,y,P),[0 Tprad],y1,opts); % warming up
    y1=Y1(end,:);
    Tpwarm=(0:h:Tprad);
    [~,Y1] = ode45(@(t,y)sistema(t,y,P),Tpwarm,y1,opts); % warming up
    y1=Y1(end,:);
    y0P(np,:)=y1;    
    toc;
end

RR=0;
RY=0;
% DR_row=[];
DR_col=[];
%Dtk=[];
Dtk=zeros(nP,L-Ndel*Ntau-1);
Rk=[];
% Dtv=[];

 for np=1:nP
     tic;
     np_v=np
     %tic;
    Iv=Ivt(1,np); % setting parameter value for learning;
y1=y0P(np,:);
[X0,Y0,Z0]=HMR_signal_v2(h,L,y1,Iv);
Zt=am*(X0-Zmn)-1; % renormalized signal y(t)
szZ=size(Zt)
%parm=Iz(1,np);
parm=Iv/Imx; % renormalized parameter;
XX=zeros(L-Ndel*Ntau-1,Ndel+1); % Mano
LL=size(XX,1)
%XX=Z(1,Ndel*Ntau+1:L-1).'; % y(t) % Kesto
% Constructing embedding delayed space:
% [y(t-tau),y(t-2*tau),..,y(t-Ndel*tau)]
for nd=1:Ndel+1 % Mano (1:Ndel -> Kesto)
    XX(:,nd)=Zt(1,(Ndel-(nd-1))*Ntau+1:L-(nd-1)*Ntau-1).'; % Mano
end
% including the parameter b:
XX=[XX, parm*ones(LL,1)]; % size(XX,2) = m+1 = Ndel+2;
Dt=(Zt(1,Ndel*Ntau+2:L)-Zt(1,Ndel*Ntau+1:L-1))/h; % dynamics of time derivative dy/dt
    %Dtk=[Dtk; Dt];
    Dtk(np,:)=Dt;
%    Dtv=[Dtv Dt];

R=ones(LL,1); % Kesto
for nd=1:degree
    R=[R,generateChebyshevs(XX,nd)];  
end
R=R.'; % This matrix will be used for optimization;
LengthW=size(R,1)
Lr=size(R,2);

Rk=[Rk R];

DR=zeros(LengthW,1);
% DR(1,1)=0.0 % since derivative of 1 is zero; (!!!)
sumt=1; % Kesto
for nd=1:degree
    DRnd=generateChebyshevsDiffLyap(XX,nd); 
    ld=length(DRnd);
    DR(sumt+1:sumt+ld,1)=DRnd;
    sumt=sumt+ld;
end
% Composing the matrices for optimization problem:
DR_col=[DR_col DR];
%DR_row=[DR_row; DR.'];
RR=RR+R*R.'; 
szDt=size(Dt)
szR=size(R)

RY=RY+R*Dt.'; % r.h.s. of optimization problem
   toc;
  end % for np=1:nP

A=[RR DR_col;
   DR_col.' zeros(nP,nP)]; % Mano;
%A=RR; % without Lc
 d=size(A,1)
 A=A+bet*eye(d);
B1=[RY; lam(1,:).'];
%B1=RY; % without Lc
W1=(A\B1).';
disp('Lagrange multipliers')
%W=W1; % without Lc
disp(W1(1,end-nP+1:end))
W=W1(1,1:end-nP);
%-----------------------------------------------------------
% Computing errors of learning:
Error=zeros(nP,1);
for np=1:nP
     Dtn=Dtk(np,1:Lr);
     Rn=Rk(:,1+(np-1)*Lr:np*Lr);
     Error(np,1)=sqrt(mean((Dtn-W*Rn).^2))/sqrt(mean(Dtn.^2)); 
end
%-----------------------------------------------------------
disp('Fitting Errors:')
disp(Error)
disp('Conditional Lyapunov Exponents:')
disp(W*DR_col)

%ur=pr
%% Computing bifurcation diagram
tpr=2000; % replication interval
Tpr=(0:h:tpr);
Npr=length(Tpr); % length of replication array
%Tpr0=1000; % initial warming interval
Tpr0=1000; % initial warming interval
Npr0=round(Tpr0/h);% number of steps in initial warming interval
r=0.006;
nu=4;
kap=-1.6;

P.r=r; 
P.nu=nu;
P.kap=kap;

opts = odeset('RelTol',1e-7,'AbsTol',1e-10);

disp("Warming up:");
tic;
Tprad=1000;
XPMaxArray=[]; % array of max of x(t) vs par-ter I;
XP0=zeros(1,Ndel*Ntau+1); % inermediate initial conditions for NGRC DDE
countv=0;
    for nb=1:Nb % scanning bifurcation parameter I;
        countv=countv+1;
    if countv==1
        toc;
        countv=0;
        nb_v=nb    
        tic;
    end
    %Iv=pM(1,Nb-nb+1);  
    Iv=pM(1,nb);  
    parm=Iv/Imx; % renormalized parameter I;   
    XP=zeros(1,Npr+Npr0); % array of solution of NGRC DDE (will be computed below)
    %----------------------------------------------------------------
 if nb==1
    %bv=pbv(1,np); % setting parameter value for learning;
    P.I=Iv;
    [~,Y1] = ode45(@(t,y)sistema(t,y,P),[0 Tprad],y0,opts); % warming up
    y1=Y1(end,:);
    Tpwarm=(0:h:Tprad);
    [~,Y1] = ode45(@(t,y)sistema(t,y,P),Tpwarm,y1,opts); % warming up
    y1=Y1(end,:);    
    % Initial conditions for warming:
    X0P=Y1(end-Ndel*Ntau:end,1).'; % X(t)
    X0P=am*(X0P-Zmn)-1; % renormalizations of x(t) for initial conditions for NGRC DDE
    %----------------------------------------------------------------
    XP(1,1:Ndel*Ntau+1)=X0P; % initial conditions; (taken from exact original system);
 else
    XP(1,1:Ndel*Ntau+1)=XP0(1,:); % intermediate initial conditions;
 end % if nb==1
generate_all_exp_matreces(m+1,degree); % generating the exponents matrices; % Mano

load("matr.mat","matr");
xx=zeros(1,m+1); 
% Integration of the NGRC DDE system:
for j=Ndel*Ntau+1:Npr-1+Npr0    
    xx(1,1:Ndel+1)=XP(1,j-(0:Ndel)*Ntau);
    xx(1,m+1)=parm;
    r=ones(LengthW,1);
    sumt=1;
    for nd=1:degree
    matrnd=matr(nd).f;
    rn=generateCheb_feat_vect(xx,matrnd,nd); % Kesto   
    ln=length(rn);
    r(sumt+1:sumt+ln,1)=rn;
    sumt=sumt+ln;
    end
    XP(1,j+1)=XP(1,j)+h*(W*r);
    % Restriction:
    % if(abs(XP(1,j+1))>1)
    % XP(1,j+1)=sign(XP(1,j));       
    % end
end
% intermediate final state of the system:
XP0=XP(1,Npr0-Ndel*Ntau:Npr0); % will be used as intermediate initial conditions
XP=zeros(1,Npr+Npr0); % array of solution of NGRC DDE (will be computed below)
TP=(0:Npr+Npr0-1); % time (-index) array;
XP(1,1:Ndel*Ntau+1)=XP0(1,:); % intermediate initial conditions;
xx=zeros(1,m+1); % Mano

% Integration of the NGRC DDE system:
for j=Ndel*Ntau+1:Npr-1+Npr0    
    xx(1,1:Ndel+1)=XP(1,j-(0:Ndel)*Ntau); % delayed variables;
    xx(1,m+1)=parm; % parameter;
    r=ones(LengthW,1);
    sumt=1;
    for nd=1:degree
    matrnd=matr(nd).f;
    rn=generateCheb_feat_vect(xx,matrnd,nd); % Kesto   
    ln=length(rn);
    r(sumt+1:sumt+ln,1)=rn; % forming the feature vector;
    sumt=sumt+ln;
    end
    XP(1,j+1)=XP(1,j)+h*(W*r); % integration;
    % Restriction:
    % if(abs(XP(1,j+1))>1)
    % XP(1,j+1)=sign(XP(1,j));       
    % end
end
% Extracting time series for finding the maxima:
XPr=XP(1,Npr0+1:Npr0+Npr); % the warming interval is omitted;
TPr=TP(1,Npr0+1:Npr0+Npr);
jXP=islocalmax(XPr(1,:)); % finding the peaks;
jXP=find(jXP); % indices of maxima
%if countv==9
   sz_jXP=size(jXP)
%end

nj=size(jXP,2);
% Filling the biurcation array:
    if nj>0
        %Iv=pM(1,Nb-nb+1); % the parameter is not normalized here; 
        Iv=pM(1,nb); 
        XPj=XPr(jXP).';
        TPj=TPr(jXP).';
        %dTPj=h*(TPj(2:end,1)-TPj(1:end-1,1));
        dTPj=TPj(2:end,1)-TPj(1:end-1,1); % time intervals between peaks;
        XPj=XPj(1:end-1,1);
        XPj=(XPj+1)/am+Zmn; % values of maxima in original normalization 
        Ivj=XPj*0+Iv;      
        %XPMaxArray=[[Ivj XPj dTPj]; XPMaxArray];
        XPMaxArray=[XPMaxArray; [Ivj XPj dTPj]]; % variables for bifurcation diagram;
    end
       
    end % for nb=1:Nb
   
toc;

%save('result_ngrc_witout_Lc_v6.mat','XPMaxArray');
save('result_ngrc_v6t.mat','XPMaxArray'); % bet=1e-6;
%save('result_ngrc_v7t.mat','XPMaxArray'); % bet=1e-7;

plot_bif_tot_v; % plot the result;


function dy = sistema(~,y,P)
% HMR neuron system:
dy=[y(2)-y(1)^3+3*y(1)^2-y(3)+P.I;
    1-5*y(1)^2-y(2);
    P.r*(P.nu*(y(1)-P.kap)-y(3))];
end
