%-------------------------------------------------------------
%--------------------------------------------------------------------
L_ex=0.91; % exact LE of Lorenz system
% Maximal period of integration:
% tpr=25000;
%-------------------------------------------------------------
%--------------------------------------------------------------------
% Parameters:
% h=0.01; tau=0.15; bet=1e-7; degree=7; Ndel=5;
% NTpr=500; (number of experiments);
% Lc_tot=[(-1:-1:-12),[-50,-100]]; % range of Lc;
%-------------------------------------------------------------
% Statistics of escape time from chaotic attractor;
%---------------------------------------------------------------------------
% Results without Lc, with standard NDRC;
 load('Tesc_without_Lc_N500.mat'); % Tpr=25000;
     Tesc=Tesc_tot;
     Lct=Lc_tot;
% Results with Lc included:    
    load('Tesc_attr_Lc_tot_N500.mat'); % Tpr=25000;
 
    Tesc=[Tesc; Tesc_tot];
    Lct=[Lct Lc_tot];

 Tesc_tot=Tesc;
 Lc_tot=Lct;
%--------------------------------------------------------------------
NTpr=size(Tesc_tot,2)
NLc=size(Tesc_tot,1)
for nTpr=1:NTpr
    Tesc_tot(:,nTpr)=flip(Tesc_tot(:,nTpr));
end
Lc_tot=flip(Lc_tot(1,:));
Lct=Lc_tot;

Tesct=[Lct; Tesc_tot.'];
Lct=Lct(1,1:end-1);
szL=size(Lct,2)
%-----------------------------------------------------------------
% Horizontal line for maximal duration of integration:
Hv=[(1:szL+1).', 2.5e4*ones(szL+1,1)]; 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
Tesc_tot_Lc=Tesct;
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Short-term valid prediction statistics; without Lc;
load('Tvp_Lc_N500.mat'); % % Ndel=5; bet=1e-7; degree=7; 

Tvpt=Tvp;
Lct=Lc_tot;

NTpr=size(Tvp,2)
NLc=size(Tvp,1)
for nTpr=1:NTpr
    Tvp(:,nTpr)=flip(Tvp(:,nTpr));
end
Lc_tot=flip(Lc_tot(1,:));
szLc=size(Lc_tot)
szTvp=size(Tvp)

Tvpt=[Lc_tot; L_ex*Tvp.'];
sz_Lc=size(Lct)
%-----------------------------------------------------------------
Lct=Tvpt(1,:);
%--------------------------------------------------------------------
% Without Lc:
load('Tvp_without_Lc_N500.mat'); % Ndel=5; bet=1e-7; degree=7;
Tvp=L_ex*Tvp.';
%--------------------------------------------------------------------
Tvpt=[Tvpt, [0; Tvp]];
sz_Tvp=size(Tvp)
%------------------------------------------------------------------------------
%-----------------------------------------------------------------
Lct=Tvpt(1,1:end-1);
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
wd=9;
hd=9;

flab=10; % Font size of text on the graph;

figure
% Plotting statistics of valid prediction time;
subplot1=subplot(2,1,1);
%boxchart(Tvpt(2:end,:),'MarkerStyle','none'); 
boxchart(Tvpt(1:end,:)); % valid prediction time statistics vs Lc;
hold(subplot1,'on');
set(gca,'TickLabelInterpreter','latex');
set(gca,'XTickLabel',{Lct,'SA'});
ylim([0,17]);
ylabel('$\Lambda T_\mathrm{vp}$','fontsize',flab,'interpreter','latex');
text(0.5,15,'(a)','FontSize',flab);
hold(subplot1,'off');

% Plotting statistics of escape time from chaotic attractor;
subplot2=subplot(2,1,2);
hold(subplot2,'on');
%boxchart(Tesct(2:end,:),'MarkerStyle','none'); % probversion
boxchart(Tesct(2:end,:)); % escape time statistics vs Lc;
plot(Hv(:,1),Hv(:,2),'r-','LineWidth',0.1); % highest limit line
set(gca,'TickLabelInterpreter','latex');
set(gca,'XTickLabel',{Lct,'SA'});
ylim([1 50000]);
ax = gca;
ax.YAxis.Scale ="log";
ytick=[1 100 1e4];
yticklab={'10^{0}','10^{2}','10^{4}'};
set(subplot2,'Ytick',ytick,'YTickLabel',yticklab,'TickLabelInterpreter','tex');
text(0.5,10000,'(b)','FontSize',flab);
xt=5.5;
yt=5.0;
% At Lc=[-8,-9] there is no escape;
text(xt,yt,'$\{$','Rotation',270,'FontSize',20,'Interpreter','latex');
text(xt,yt+0.5,'no escape','Rotation',90,'FontSize',10,'Interpreter','latex');

xlabel('$\lambda_{\perp}$','Interpreter','latex'); % Lc
ylabel('$T_\mathrm{esc}$','fontsize',flab,'interpreter','latex');

ax2=gca;
v2=ax2.Children
ax2.Children=ax2.Children([1 2 3 5 4]);
hold(subplot2,'off');

set(gcf,'Units','centimeters');
set(gcf, 'PaperSize', [wd,hd]);
set(gcf,'Position',[2,2,wd,hd]);
savefig('Lor_Tesc_Tvp_Lc.fig')
saveas(gcf,'Lor_Tesc_Tvp_Lc','epsc')