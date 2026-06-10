

% Result for NGRC DDE system learned from HMR system:
load('result_ngrc_v6t.mat','XPMaxArray'); % without restriction;

h=0.05; % step of integration
dn=1;
bv=XPMaxArray(1:dn:end,1); % array of parameters;
Ymx=h*XPMaxArray(1:dn:end,3); % array of time intervals;

nP=5;
% Set of system parameters for which the learning was performed;
pbv=linspace(3.0,3.5,nP); % with nP=5; 
Yv=[0; 100];

% Result for exact HMR system:
load('result_exact_v.mat','XPMaxArray');

bvt=XPMaxArray(:,1); % array of parameters;
Ymxt=XPMaxArray(:,3); % array of time intervals;

fnum=10;
flab=10;

wd=9;
hd=9;

xa=3.01;
ya=93;

figure
subplot1=subplot(2,1,1)
hold on
% Plotting bifurcation diagram for original HMR system;
plot(bvt,Ymxt,'b.','MarkerSize',1); % exat HMR
%plot(bv,Ymx,'.','MarkerSize',2,'Color',[1 0.5 0]); % NGRC from HMR
% Plotting vertical lines of parameters for which the 
% learning was performed;
for np=1:nP
plot([pbv(1,np); pbv(1,np)],Yv(:,1),'r--','LineWidth',1);
end
% a = get(subplot1,'XTickLabel');  
% a1={'3','3.1','3.2','3.3','3.4','3.5'};
% set(subplot1,'XTickLabel',a1,'fontsize',fnum)
xlim([3.0 3.5]);
ylim([0 100]);
ylabel('$\Delta T_{\mathrm{n}}$','fontsize',flab,'interpreter','latex');
ax1 = gca;
ax1.XAxis.FontSize = flab;
ax1.YAxis.FontSize = flab;
text(xa+0.39,ya,'original','FontSize',flab);
text(xa,ya,'(a)','FontSize',flab);
hold off
subplot2=subplot(2,1,2)
hold on
% Plotting bifurcation diagram for NGRC DDE system
% that was learned from the original HMR system;
plot(bv,Ymx,'*','MarkerSize',1,'Color',[1 0.5 0]); % NGRC from HMR
for np=1:nP
plot([pbv(1,np); pbv(1,np)],Yv(:,1),'r--','LineWidth',1);
end
% a = get(subplot2,'XTickLabel')
% a2={'3','3.1','3.2','3.3','3.4','3.5'};
% set(subplot2,'XTickLabel',a2,'fontsize',fnum)
ax2 = gca;
ax2.XAxis.FontSize = flab;
ax2.YAxis.FontSize = flab;
text(xa+0.38,ya,'reservoir','FontSize',flab);
text(xa,ya,'(b)','FontSize',flab);
ylabel('$\Delta T_{\mathrm{n}}$','fontsize',flab,'interpreter','latex');
xlabel('$I_{\mathrm{ext}}$','fontsize',flab,'interpreter','latex');
xlim([3.0 3.5]);
hold off

set(gcf,'Units','centimeters');
set(gcf, 'PaperSize', [wd,hd]);
set(gcf,'Position',[2,2,wd,hd]);
savefig('HMR_bif.fig')
saveas(gcf,'HMR_bif','epsc')

