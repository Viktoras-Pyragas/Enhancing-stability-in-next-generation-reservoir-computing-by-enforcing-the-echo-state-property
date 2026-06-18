%----------------------------------------------------------------
% main_Fig_MG.m
%-------------------------------------------------------------------------------------------
% This file plots the Figure for the Mackey-Glass system;
% The script compares the original Mackey-Glass trajectory with the
% autonomous NGRC/reservoir trajectory stored in Xtot_record.mat.
%-------------------------------------------------------------------------------------------
close all
format long;
h=0.01;
lw2=0.5;
%--------------------------------------------------
% Basic plotting parameters used throughout the figure.
color_ex=[0 0.3 1]; % color of original solution
color_rc=[1 0.5 0]; % color of solutions of NGRC
segLength = 80; % period of the same color tone
% Color tones parameters (a1,a2):
% age = k/N;
% alpha = a1 + a2*age;   % older/fainter to newer/darker
% original values:
% a1=0.03;
% a2=0.25;
% probe values:
a1=0.03;
a2=0.5;
lwd=0.7; % linewidth for phase portraits;
%--------------------------------------------------
% Dynamics of the exact system and of NGRC after learning;
%--------------------------------------------------
% Xtot=[Tpr1 U1 XP1];
% save('XPt_record2.mat','Xtot');
% Tpr1 -> time sampling;
% U1 -> exact time series;
% XP1 -> time series of NGRC;
%--------------------------------------------------
% epsnorm=0.0; Lc_tot=-1.0; h=0.01; tau=0.5; TL=1000;
% bet=1e-7; degree=6; Ndel=4; Tprad_tot=0;
% Load the matrix Xtot, whose columns are time, original solution, and
% NGRC/reservoir solution.
load('Xtot_record.mat');
%--------------------------------------------------
%Xtot=[Tpr1 U1 XP1]; % array from the mat-file
  Tpr1=Xtot(:,1); % time sampling
  Tpr1=Tpr1-Tpr1(1,1); % time sampling starts from zero
  U1=Xtot(:,2); % dynamics of original system
  XP1=Xtot(:,3); % dynamics of autonomous NGRC system
  
%-----------------------------------------------------------    
% Create the final multi-panel figure.
figure

% Coordinates for panel labels (a)-(e).
X1a=[0.5 0.9];
Xat=X1a(1);
Yat=X1a(2);

X1b=[0.7 0.9];
Xbt=X1b(1);
Ybt=X1b(2);

X1c=[0.7 0.9];
Xct=X1c(1);
Yct=X1c(2);

X1d=[-0.9 1.3];
Xdt=X1d(1);
Ydt=X1d(2);

X1e=[1.6 10];
Xet=X1e(1);
Yet=X1e(2);

% Panel (a): time-domain comparison of the original and NGRC signals.
subplot(3,2,[1 2]);
wideAx=gca;
axes(wideAx);
hold on
% solution of original system:
plot(Tpr1,U1,'-','color',color_ex,'LineWidth',lw2) 
% solution of NGRC system:
plot(Tpr1,XP1,'-','color',color_rc,'LineWidth',lw2) 
hold off
xlim([0,100]);
ylim([-1.1 1.1]);
xlabel('time','Interpreter','latex','FontSize',10);
ylabel('$u$, $\nu$','Interpreter','latex','FontSize',10);
text(Xat,Yat,'(a)','FontSize',10);
leg1=legend('original', 'reservoir', 'Orientation', 'horizontal');
leg1.Position = [0.49 0.95 0.1 0.01];
set(wideAx,'XScale','linear');
set(wideAx,'YScale','linear');

set(wideAx,'XTick',[0 20 40 60 80 100]);
set(wideAx,'YTick',[-1 0 1]) 
set(wideAx,'XTickLabel',{'0','20','40','60','80','100','Interpreter','latex','FontSize',10})
set(wideAx,'YTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
box on
%------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
% Compute delayed-coordinate data for the phase portraits.
% Xnow_ex and Xdel_ex belong to the original system.
% Xnow and Xdel belong to the NGRC/reservoir system.
[Xnow,Xdel,Xnow_ex,Xdel_ex] = plot_MG_dyn(Xtot);
lw2=0.1;
% Panel (b): delayed phase portrait of the original Mackey-Glass system.
subplot(3,2,3);
ha=gca;
axes(ha);
hold on;
%-------------------------------------------------------------
N=size(Xnow_ex,1);
%segLength = 80;
x1=Xnow_ex;
x2=Xdel_ex;

% Draw the phase trajectory in short segments so older and newer parts can
% be shown with different transparency levels.
for k = 1:segLength:N-segLength
    age = k/N;
    alpha = a1 + a2*age;   % older/fainter to newer/darker
      
    plot(x1(k:k+segLength), x2(k:k+segLength), ...
        'Color',[color_ex alpha], ...
        'LineWidth',lwd);
end
hold off
%-------------------------------------------------------------
set(ha,'XTickLabel','');
set(ha,'XScale','linear');
set(ha,'YScale','linear');
xlabel('$u(t)$','Interpreter','latex','FontSize',10);
ylabel('$u(t-\tau)$','Interpreter','latex','FontSize',10);
text(Xbt,Ybt,'(b)','FontSize',10);

set(ha,'XTick',[-1 0 1]);
set(ha,'YTick',[-1 0 1]) 
set(ha,'XTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})
set(ha,'YTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
xlim([-1.1,1.1])
ylim([-1.1,1.1])
 box on

% Panel (c): delayed phase portrait of the autonomous NGRC/reservoir system.
subplot(3,2,4);
ha=gca;
axes(ha);
hold on;
%-------------------------------------------------------------
N=size(Xnow,1)
%segLength = 80;
x1=Xnow;
x2=Xdel;

% Use the same segmented transparency style as in panel (b), but with the
% NGRC/reservoir color.
for k = 1:segLength:N-segLength
    age = k/N;
    alpha = a1 + a2*age;   % older/fainter to newer/darker
        
    plot(x1(k:k+segLength), x2(k:k+segLength), ...
        'Color',[color_rc alpha], ...
        'LineWidth',lwd);
end
hold off
%-------------------------------------------------------------
xlabel('$\nu(t)$','Interpreter','latex','FontSize',10);
ylabel('$\nu(t-\tau)$','Interpreter','latex','FontSize',10);
text(Xct,Yct,'(c)','FontSize',10);

set(ha,'XTick',[-1 0 1]);
set(ha,'YTick',[-1 0 1]) 
set(ha,'XTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})
set(ha,'YTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
xlim([-1.1,1.1])
ylim([-1.1,1.1])
 box on
        % Compute probability density functions for both signals.
        [bin_centers,pdf_ngrc,pdf_exact] = plot_MG_pdf(Xtot);
% Panel (d): probability density comparison.
subplot(3,2,5);
ha=gca;
axes(ha);
plot(bin_centers, pdf_exact,'color',color_ex, 'LineWidth', lw2);
hold on;
plot(bin_centers, pdf_ngrc,'color', color_rc, 'LineWidth', lw2);
hold off;

set(ha,'XTickLabel','');
set(ha,'XScale','linear');
set(ha,'YScale','linear');
xlabel('$u,\nu$','Interpreter','latex','FontSize',10);
ylabel('PDF','Interpreter','latex','FontSize',10);
text(Xdt,Ydt,'(d)','FontSize',10);

set(ha,'XTick',[-1 0 1]);
set(ha,'YTick',[0 0.5 1 1.5]) 
set(ha,'XTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})
set(ha,'YTickLabel',{'0','0.5','1.0','1.5','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
xlim([-1,1])
ylim([0,1.5])
 box on
%------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
% Compute the power spectral density of the original and NGRC signals.
[PSD_exact,PSD_ngrc,f] = plot_MG_psd(Xtot);
% save('MG_psd.mat','f','PSD_exact','PSD_ngrc');
%load('MG_psd.mat');

% Panel (e): frequency-domain comparison using logarithmic power scale.
subplot(3,2,6);
ha=gca;
semilogy(f, PSD_exact,'color',color_ex, 'LineWidth', lw2);
hold on;
semilogy(f, PSD_ngrc,'color',color_rc, 'LineWidth', lw2);
hold off;

ylabel('Power','Interpreter','latex','FontSize',10);
xlim([0,2])
%xlim([0,5])
%ylim([1e-8,1])
ylim([1e-7,100])
text(Xet,Yet,'(e)','FontSize',10);

set(ha,'XTick',[0 1 2]);
set(ha,'XTickLabel',{'0','1','2','Interpreter','latex','FontSize',10})
xlabel('Frequency','Interpreter','latex','FontSize',10);
set(ha,'XScale','linear');
set(ha,'YScale','log');

set(ha,'YTick',[1e-7 1e-3 10]) 
set(ha,'YTickLabel',{'10^{-7}','10^{-3}','10^{1}','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
box on
%------------------------------------------------------------------------------------------------
% Set the physical figure size and save the final result as an EPS file.
set(gcf,'Units','centimeters')
set(gcf, 'PaperSize', [8.5 11]);
set(gcf,'Position',[2,2,8.5,11]);
saveas(gcf,'Fig_MG','epsc')
