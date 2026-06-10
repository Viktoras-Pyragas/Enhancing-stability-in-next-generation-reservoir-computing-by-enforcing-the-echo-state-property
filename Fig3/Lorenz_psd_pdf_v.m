close all
clearvars

format long;

load('Xvt4t.mat','Ut'); % exact Lorenz system; h=0.01; NT=1e7;

load('Tesc_XP5t.mat','Xt'); % NGRC DDE Lorenz system; % h=0.01; NT=1e7;

Tv=Ut(:,1);
x=Ut(:,2);

M=19.3667;  % max computed from a very long (T=10000) tme interval;
x=x/M;

figure
hold on
plot(Ut(:,1),Ut(:,2)/M,'b-');
plot(Xt(:,1),Xt(:,2),'r--');
hold off


Tv=Tv-Tv(1,1);
Nt=size(Tv,1)
%ur=pr
Kt=1000; % number of segments
dt=Tv(2,1)-Tv(1,1) % time step
fs=1/dt

% Spectra from exact signal:
[Power1, Frequency1] = WelchPowerSpectralDensity(x, [], Hann(round(Nt/Kt)), 0.0, fs);

Tv=Xt(:,1);
x=Xt(:,2);
%x=x/M;

Tv=Tv-Tv(1,1);
Nt=size(Tv,1)
%ur=pr
dt=Tv(2,1)-Tv(1,1) % time step
fs=1/dt

[Power2, Frequency2] = WelchPowerSpectralDensity(x, [], Hann(round(Nt/Kt)), 0.0, fs);

%--------------------------------------------------------------------
% Computing pdf for exact Lorenz model;
%--------------------------------------------------------------------
hT=0.01;

load('Xvt4t.mat','Ut'); % time series of exact Lorenz system;
X=Ut(:,2);

Ut=sort(X(:,1),'ascend');
Ut=Ut/M;

Nx=size(Ut,1)
%ur=pr
dx=0.02;
xpt=(-1:dx:1);
szx=size(xpt,2)
count_tot=zeros(szx-1,2);
for nx=1:szx-1
   edges=[xpt(1,nx) xpt(nx+1)];
   [counts,edges] = histcounts(Xt,edges);
   count_tot(nx,1:2)=[xpt(1,nx) counts];  
end
Nrm=dx*sum(count_tot(:,2))
count_tot(:,2)=count_tot(:,2)/Nrm;

count_ex=count_tot;

%--------------------------------------------------------------------
% Computing pdf for ngrc dde learned from Lorenz model;
%--------------------------------------------------------------------
load('Tesc_XP5t.mat','Xt'); % time series of NGRC DDE system learned from Lorenz system;
Xt=sort(Xt(:,2),'ascend');

Nx=size(Xt,1)
dx=0.02;
xpt=(-1:dx:1);
szx=size(xpt,2)
count_tot=zeros(szx-1,2);
for nx=1:szx-1
   edges=[xpt(1,nx) xpt(nx+1)];
   [counts,edges] = histcounts(Xt,edges);
   count_tot(nx,1:2)=[xpt(1,nx) counts];  
end
Nrm=dx*sum(count_tot(:,2))
count_tot(:,2)=count_tot(:,2)/Nrm;
%--------------------------------------------------------------------
%--------------------------------------------------------------------
wd=9;
hd=9;
fnt=10;
figure
subplot1=subplot(2,1,1);
semilogy(Frequency1,abs(Power1),'color', 'b','LineStyle', '-','LineWidth',1); % exact Lorenz
hold on
semilogy(Frequency2,abs(Power2),'color', [1 0.5 0],'LineStyle', '-','LineWidth',1); % HMR from NGRC DDE
%hold on
setx=0.2;
sety=0.7e-2;
text(setx,sety,'(a)','Color','black','FontSize',fnt);
lgd1=legend({'original','reservoir'});
legend boxoff
%posv=[0.719880950575783   0.832764391935278   0.166785716090884   0.075476192020235];
posv=[0.67   0.832764391935278   0.166785716090884   0.075476192020235];
lgd1.Position=posv;

xlim([0 15]);
xlabel('Frequency',Interpreter='latex');
ylabel('Power',Interpreter='latex');
hold off
subplot2=subplot(2,1,2);
plot(count_ex(:,1),count_ex(:,2),'color', 'b','LineStyle', '-','LineWidth',1);
hold on
plot(count_tot(:,1),count_tot(:,2),'color', [1 0.5 0],'LineStyle', '-','LineWidth',1);
setx=-0.97;
sety=0.88;
text(setx,sety,'(b)','Color','black','FontSize',fnt);
lgd2=legend({'original','reservoir'});
posv=[0.67   0.358478677649564   0.166785716090884   0.075476192020235];
lgd2.Position=posv;
legend boxoff
xlabel('$\nu$, $u$',Interpreter='latex');
ylabel('PDF',Interpreter='latex');
hold off
set(gcf,'Units','centimeters');
set(gcf, 'PaperSize', [wd,hd]);
set(gcf,'Position',[2,2,wd,hd]);
savefig("Lor_psd_pdf.fig")
saveas(gcf,'Lor_psd_pdf','epsc')
%--------------------------------------------------------------------
%--------------------------------------------------------------------
