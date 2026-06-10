%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Tn-Tn1 map:
%------------------------------------------------------------------------
% Xt=[Tpr.' XP];
% save('XP_dyn.mat','Xt');
%----------------------------------------------------------
%load('XP_dyn.mat'); % time series from NGRC DDE learned by HMR system
load('XP_dyn_v2t.mat'); % time series from NGRC DDE learned by HMR system

% Computing the map of Tn1=f(Tn)
jmaxt=find(islocalmax(Xt(:,2)));
XF=[Xt(jmaxt,1) Xt(jmaxt,2)]; % series of maxima (Tn, Xn)

XF=[XF(2:end,1)-XF(1:end-1,1) XF(1:end-1,2)]; % series of  (dTn, Xn)

TXFn=[XF(1:end-1,1) XF(2:end,1)]; % Map series of (dTn, dTn1)
%----------------------------------------------------------
% Ut=[Tpr.' U];
% save('Ut_dyn.mat','Ut');
%load('Ut_dyn.mat'); % time series from exact HMR system
load('Ut_dyn_v2t.mat'); % time series from exact HMR system

% Renormalizing time series;
Zmx=max(Ut(:,2));
Zmn=min(Ut(:,2));
am=2/(Zmx-Zmn);

Ut(:,2)=am*(Ut(:,2)-Zmn)-1; % Renormalized signal, (-1<=Ut<=1);

% Computing the map of Tn1=f(Tn)
jmaxt=find(islocalmax(Ut(:,2)));
UF=[Ut(jmaxt,1) Ut(jmaxt,2)]; % series of maxima (Tn, Un)

UF=[UF(2:end,1)-UF(1:end-1,1) UF(1:end-1,2)]; % series of (dTn, Un)

TUFn=[UF(1:end-1,1) UF(2:end,1)];  % Map series of (dTn, dTn1)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% pdf for Tn
%------------------------------------------------------------------------

% Xt=[Tpr.' XP];
% save('XP_dyn.mat','Xt');

Ndx=200; % number of time segments for PDF distribution;
%----------------------------------------------------------
load('XP_dyn_v2t.mat'); % time series from NGRC DDE learned by HMR system

% Computing the local maxima of time series
jmaxt=find(islocalmax(Xt(:,2)));
XF=[Xt(jmaxt,1) Xt(jmaxt,2)];

% Series of time intervals between spikes (dTn=Tn1-Tn):
Tn_ngrc=XF(2:end,1)-XF(1:end-1,1); 

%----------------------------------------------------------
% Ut=[Tpr.' U];
% save('Ut_dyn.mat','Ut');
load('Ut_dyn_v2t.mat'); % time series from exact HMR system

% Renormalizing time series;
Zmx=max(Ut(:,2));
Zmn=min(Ut(:,2));
am=2/(Zmx-Zmn);

Ut(:,2)=am*(Ut(:,2)-Zmn)-1; % Renormalized signal, (-1<=Ut<=1);

% Compuiting the local maxima of time series
jmaxt=find(islocalmax(Ut(:,2)));
UF=[Ut(jmaxt,1) Ut(jmaxt,2)];

% Series of time intervals between spikes (dTn=Tn1-Tn):
Tn_ex=UF(2:end,1)-UF(1:end-1,1);
%--------------------------------------------------------------------
% Computing pdf for exact HMR model;
%--------------------------------------------------------------------
X=Tn_ex;
% Evaluating the range of dTn 
xmn=min(X);
xmx=max(X);
dxt=xmx-xmn;

% xmn1=xmn;
% xmx1=xmx;
% dxt1=dxt;

% Setting the range parameters
xmn1=9.5;
xmx1=76.5;
dxt1=xmx1-xmn1;

Xt=sort(X(:,1),'ascend');

Nx=size(Xt,1)
dx=dxt1/Ndx
xpt=(xmn1:dx:xmx1);
szx=size(xpt,2)
count_tot=zeros(szx-1,2);
% Computing the PDF:
for nx=1:szx-1
   edges=[xpt(1,nx) xpt(nx+1)];
   [counts,edges] = histcounts(Xt,edges);
   count_tot(nx,1:2)=[xpt(1,nx) counts];  
end
% Normalization:
Nrm=dx*sum(count_tot(:,2))
count_tot(:,2)=count_tot(:,2)/Nrm;

count_ex=count_tot;

%--------------------------------------------------------------------
% Computing pdf for ngrc dde learned from HMR model;
%--------------------------------------------------------------------
Xt=Tn_ngrc;
Xt=sort(Xt(:,1),'ascend');
% Evaluating the range of dTn 
xmn=min(Xt)
xmx=max(Xt)
dxt=xmx-xmn

% xmn1=xmn;
% xmx1=xmx;
% dxt1=dxt;

% Setting the range parameters
xmn1=9.5;
xmx1=76.5;
dxt1=xmx1-xmn1;

Nx=size(Xt,1)
dx=dxt1/Ndx;
xpt=(xmn1:dx:xmx1); % the range of dTn segments
szx=size(xpt,2)
count_tot=zeros(szx-1,2);
% Computing the PDF:
for nx=1:szx-1
   edges=[xpt(1,nx) xpt(nx+1)];
   [counts,edges] = histcounts(Xt,edges);
   count_tot(nx,1:2)=[xpt(1,nx) counts];  
end
% Normalization:
Nrm=dx*sum(count_tot(:,2))
count_tot(:,2)=count_tot(:,2)/Nrm; 
%--------------------------------------------------------------------
%--------------------------------------------------------------------
% Plotting the general graph;
%--------------------------------------------------------------------
%--------------------------------------------------------------------
wd=9;
hd=9;
fnt=10;
figure
subplot1=subplot(2,1,1);
p1=plot(TUFn(:,1),TUFn(:,2)); % exact HMR, (dTn,dTn1)
hold on
p2=plot(TXFn(:,1),TXFn(:,2)); % HMR from NGRC DDE, (dTn,dTn1)
p1.Marker='*';
p2.Marker='*';
p1.MarkerSize=6;
p2.MarkerSize=2;
p1.LineStyle='none';
p2.LineStyle='none';
p1.Color='b';
p2.Color=[1 0.5 0];
setx=12.0;
sety=72.0;
text(setx,sety,'(a)','Color','black','FontSize',fnt);

lgd1=legend({'original','reservoir'});
legend boxoff
posv=[0.67   0.832764391935278   0.166785716090884   0.075476192020235];
lgd1.Position=posv;

xlim([10 80]);
ylim([0 80]);
xlabel('$\Delta T_{n}$','Interpreter','latex');
ylabel('$\Delta T_{n+1}$','Interpreter','latex');
hold off
subplot2=subplot(2,1,2);
% PDF of exact HMR model
plot(count_ex(:,1),count_ex(:,2),'color', 'b','LineStyle', '-','LineWidth',1);
hold on
% PDF of NGRC DDE learned from HMR model
plot(count_tot(:,1),count_tot(:,2),'color', [1 0.5 0],'LineStyle', '-','LineWidth',1);
setx=12.0;
sety=0.18;
text(setx,sety,'(b)','Color','black','FontSize',fnt);
lgd2=legend({'original','reservoir'});
posv=[0.67   0.358478677649564   0.166785716090884   0.075476192020235];
lgd2.Position=posv;
legend boxoff
xlim([10 80]);
xlabel('$\Delta T_{n}$','Interpreter','latex');
ylabel('PDF','Interpreter','latex');
hold off
set(gcf,'Units','centimeters');
set(gcf, 'PaperSize', [wd,hd]);
set(gcf,'Position',[2,2,wd,hd]);
savefig("HMR_map_pdf.fig")
saveas(gcf,'HMR_map_pdf','epsc')
%--------------------------------------------------------------------
%--------------------------------------------------------------------
