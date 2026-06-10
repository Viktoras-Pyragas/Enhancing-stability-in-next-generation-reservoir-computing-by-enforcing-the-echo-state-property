clear all

L_ex=0.013359437831911; % exact leading LE of HMR system
% Maximal period of integration:
Tmx1e5=1e5; % tpr=1e5;
%-------------------------------------------------------------
% Parameters of computation:
% NTpr=500; beta=1e-6; h=0.05; tau=1.0; degree=8; Ndel=5;
% TL=2000; % Learning interval;
%-------------------------------------------------------------
%---------------------------------------------------------------------------

Tvp_t=[];
Lc_t=[];
load('Tvp_tot_Lc_N500_tpr1e3.mat');  
% Lc=[-0.08 -1 -2 -3 -3.2 -5 -10 -20]; bet=1e-6; tpr=1000; NTpr=500;

NTpr=size(Tvp_tot,2);
NLc=size(Tvp_tot,1);
Tvp_t=[Tvp_t; Tvp_tot];
Lc_t=[Lc_t, Lc_tot];

NLc=size(Tvp_t,1);

for nTpr=1:NTpr
    Tvp_t(:,nTpr)=flip(Tvp_t(:,nTpr));
end
Lc_t=flip(Lc_t(1,:));
Lct=Lc_t;

Tvp_tot_Lc=L_ex*Tvp_t.';
szTvpc=size(Tvp_tot_Lc)

load('Tvp_tot_without_Lc_N500_tpr1e3.mat');
% bet=1e-6; tpr=1000; NTpr=500;

szT=size(Tvp_tot)
szT_Lc=size(Tvp_tot_Lc)

Tvp_tot_Lc=[Tvp_tot_Lc, L_ex*Tvp_tot.'];
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Lc=[-0.08 -1 -2 -3 -3.2 -5 -10 -20]; 
% bet=1e-6; tpr=1e5; NTpr=1e7;
load('Tesc_tot_Lc_N500_tpr1e5.mat');
Tesc_t=Tesc_tot;
Lc_t=Lc_tot;


for nTpr=1:NTpr
    Tesc_t(:,nTpr)=flip(Tesc_t(:,nTpr));
end
Lc_t=flip(Lc_t(1,:));
Lct=Lc_t;
Lc_tot=Lct;

szL=size(Lct,2)
Hv=[(1:szL+1).', Tmx1e5*ones(szL+1,1)];
Tesc_tot_Lc=Tesc_t.';
szTvesc=size(Tesc_tot_Lc)

load('Tesc_tot_without_Lc_N500_tpr1e5.mat');
% bet=1e-6; tpr=1e5; NTpr=1e7;

szT=size(Tesc_tot)
szT_Lc=size(Tesc_tot_Lc)
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
szT=size(Tesc_tot)
szT_Lc=size(Tesc_tot_Lc)
Tesc_tot_Lc=[Tesc_tot_Lc, Tesc_tot.'];

szT=size(Tesc_tot_Lc);
n=szT(1);
m=szT(2);
for n1=1:n
    for m1=1:m
        if Tesc_tot_Lc(n1,m1)==0       
            % These solutions did not escape from chaotic attractor;
            Tesc_tot_Lc(n1,m1)=Tmx1e5;  
        end
    end
end

szTvp=size(Tvp_tot_Lc)
szTesc=size(Tesc_tot_Lc)

wd=9;
hd=9;

flab=10;
fnum=10;
rd=45;

figure
subplot1=subplot(2,1,1)
set(subplot1,'PositionConstraint','innerposition')
hold(subplot1,'on');
%boxchart(Tvp_tot_Lc,'MarkerStyle','none');
boxchart(Tvp_tot_Lc);
set(subplot1,'TickLabelInterpreter','latex');
%xtickangle(rd);
set(subplot1,'XTickLabel',{'-20','-10','-5','-3.2','-3','-2','-1','-0.08','SA'});
%xl=xlim
ylim(subplot1,[0 6]);
a = get(subplot1,'XTickLabel');  
set(subplot1,'XTickLabel',a,'fontsize',fnum)
%xlabel('$\lambda_{\perp}$','fontsize',flab,'interpreter','latex');
ylabel('$T_{\mathrm{vp}}\Lambda$','fontsize',flab,'interpreter','latex');
ax1 = gca;
ax1.XAxis.FontSize = flab;
ax1.YAxis.FontSize = flab;
hold(subplot1,'off');

subplot2=subplot(2,1,2)
set(subplot2,'PositionConstraint','innerposition')
%boxchart(Tesc_tot_Lc,'MarkerStyle','none');
boxchart(Tesc_tot_Lc);
hold(subplot2,'on');
plot(Hv(:,1),Hv(:,2),'r-','LineWidth',0.1);
set(subplot2,'TickLabelInterpreter','latex');
%xtickangle(rd);
set(subplot2,'XTickLabel',{'-20','-10','-5','-3.2','-3','-2','-1','-0.08','SA'});

ylim(subplot2,[1 1.5e5]);
a = get(subplot2,'XTickLabel');  
set(subplot2,'XTickLabel',a,'fontsize',fnum)
ytick=[1 100 1e4];
yticklab={'10^{0}','10^{2}','10^{4}'};
set(subplot2,'Ytick',ytick,'YTickLabel',yticklab,'TickLabelInterpreter','tex');
xlabel('$\lambda_{\perp}$','fontsize',flab,'interpreter','latex');
ylabel('$T_{\mathrm{esc}}$','fontsize',flab,'interpreter','latex');
% ax = gca;
% ax.YAxis.Scale ="log";

xt=6.0;
yt=20.0;
% At Lc=[-1,-2,-3] there is no escape;
%text(xt,yt,'$\{$','Rotation',270,'FontSize',50,'Interpreter','latex');
text(xt,yt,'$\Bigg\{$','Rotation',270,'FontSize',15,'Interpreter','latex');
text(xt-0.1,yt,'no escape','Rotation',90,'FontSize',9,'Interpreter','latex');

  
ax2 = gca;
ax2.XAxis.FontSize = flab;
ax2.YAxis.FontSize = flab;
ax2.YAxis.Scale ="log";
v2=ax2.Children
ax2.Children=ax2.Children([1 2 4 3]);
hold(subplot2,'off');

set(gcf,'Units','centimeters');
set(gcf, 'PaperSize', [wd,hd]);
set(gcf,'Position',[2,2,wd,hd]);
savefig('HMR_Tesc_Tvp_Lc.fig')
saveas(gcf,'HMR_Tesc_Tvp_Lc','epsc')