% Prediction of X variable of the ele QIF networksystem with a
% recunstructed model which takes into account the negative transversal
% Lyapunov exponent

% Clean the MATLAB workspace and close any figures from previous runs.
close all
clearvars

% Main model parameters.
% lam sets the desired conditional Lyapunov exponent of the reconstructed
% model. A negative value means that the reconstructed dynamics should be
% locally stable in the selected transverse direction.
%lam=-.01; % desired conditional Lyapunov exponent
lam=-1; % desired conditional Lyapunov exponent
%h=.01; % sampling time
h=.04; % sampling time
Ntau=15; % Numer of points in embedding time interval

% Load the simulated QIF network data. The script uses every second sample
% from the loaded QIF array: column 1 is time and column 2 is the signal to
% be modeled.
 load QIF_netw_ODE4_N1e5_h02_3_2.mat;
 v2=QIF(1:2:end,2);
 tm=QIF(1:2:end,1);
Lv2=length(v2);
disp('Total length of the experimental signal');
disp(Lv2); % Total number of points in the experimental signal
Lv2=length(v2);

% Choose the length of the learning interval and prediction interval, and
% define the nonlinear model size.
L=500000; % Number of ponts in the learning interval
LP=750000; % Number of points in the prediction interval
bet=2e-6; % ridge regression parameter
degree=5; % degree of nonlinearity
Ndel=9; % number of delays
m=Ndel+1; % embedding dimension

% Training signal. It is stored as a row vector because later indexing and
% feature construction use row-vector time series.
Y0=v2(1:L).';
eps=0.0; % small parameter to scale the signal in the interval (-1+eps, 1+eps) 

% Rescale the training signal to a normalized interval close to [-1, 1].
% Chebyshev polynomial bases are naturally defined on this interval.
Zmn=min(v2);
Zmx=max(v2);
am=2/(Zmx-Zmn);
Z=(1-eps)*(am*(Y0-Zmn)-1);

% Construct the delayed embedding matrix XX. The first column is the present
% value, and the following columns contain delayed values separated by Ntau
% samples.
XX=Z(Ndel*Ntau+1:L-1).';
for nd=1:Ndel
    XX=[XX, Z((Ndel-nd)*Ntau+1:L-nd*Ntau-1).'];
end
LL=size(XX,1);

% Approximate the time derivative of the normalized signal by a forward
% finite difference. This is the regression target.
D=(Z(Ndel*Ntau+2:L)-Z(Ndel*Ntau+1:L-1))/h;

% Build the full regression feature matrix using Chebyshev polynomial terms
% from degree 1 up to the selected maximum degree. The constant term is
% included first.
R=ones(LL,1);
for nd=1:degree
    R=[R,generateChebyshevs(XX,nd)];
end
R=R.';
LengthW=size(R,1);

% Build the derivative vector of the Chebyshev features. This vector is used
% to impose the desired conditional Lyapunov exponent during fitting.
DR=zeros(LengthW,1);
sm=1;
for nd=1:degree
    DRnd=generateChebyshevsDiffLyap(XX,nd);
    ld=length(DRnd);
    DR(sm+1:sm+ld)=DRnd;
    sm=sm+ld;
end

% Assemble and solve the constrained ridge-regression problem.
% RR contains the normal-equation matrix for the feature data. The extra row
% and column add a Lagrange multiplier enforcing W * DR = lam.
RR=R*R.';
A=[RR,DR];
A=[A;[DR.',0]];
d=size(A,1);
A=A+bet*eye(d);
B1=[R*D.';lam];
W1 = (A\B1).';
%W1 = (pinv(A)*B1).';

% Report fitting diagnostics: the Lagrange multiplier, relative fitting
% error, and achieved conditional Lyapunov exponent.
disp('Lagrange multiplier')
disp(W1(end))
W=W1(1:end-1);
Error=sqrt(mean((D-W*R).^2))/std(D);
disp('Fitting Error')
disp(Error)
disp('Condtional Lyapunov')
disp(W*DR)

%disp(W)
%T=linspace(0,hT*(NT-1),NT);

% Prediction after learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the prediction interval and the real normalized signal used for
% comparison with the model forecast.
Tpr=0:h:LP*h;
Npr=length(Tpr); % Prognozes masyvo ilgis
Y1=v2(L:L+Npr-1).';
XP=zeros(1,Npr);
j0=Ndel*Ntau+1; % initial point for itteration of the forrecusting model
U=am*(Y1-Zmn)-1; % original signal

% Initialize the model prediction with true data so that all required delay
% coordinates are available before autonomous prediction starts.
XP(1:Ndel*Ntau+1)=U(1:Ndel*Ntau+1);

% Precompute exponent matrices used for fast Chebyshev feature construction
% during the prediction loop.
generate_all_exp_matreces(m,degree);
load("matr.mat");
xx=zeros(1,m);
tic

% Iterate the learned differential model forward in time using Euler's
% method. At each step, the current delayed state is converted into a
% Chebyshev feature vector and multiplied by the learned weights W.
for j=j0:LP-1
    for jm=1:m
        xx(jm)=XP(j-(jm-1)*Ntau);
    end
    r=ones(LengthW,1);
    sm=1;
    for nd=1:degree
    matrnd=matr(nd).f;
    rn=generateCheb_feat_vect(xx,matrnd, nd);
    ln=length(rn);
    r(sm+1:sm+ln)=rn;
    sm=sm+ln;
    end
    XP(j+1)=XP(j)+h*(W*r);    
end
toc

% Plot the real normalized signal and the predicted signal on the same time
% axis. The red star marks the point where autonomous prediction begins.
TL=LP*h; % prediction time interval
figure
hold on
plot(TL+Tpr,U,'-b')
plot(TL+Tpr,XP,'-r')
plot(TL+Tpr(j0),XP(j0),'*r')
ylim([-1 1])

% Save time, true normalized signal, and predicted signal for later analysis.
TUV=[Tpr.', U.', XP.'];
save("TUV.mat","TUV");

