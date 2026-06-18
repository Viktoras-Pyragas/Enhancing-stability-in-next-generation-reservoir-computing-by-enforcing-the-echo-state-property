function Plot_Short_Long_Intersp_map(Qp,Qd)

[lw2,lw3,lwd,color_ex,color_rc,a1,a2,segLength,Dt,Nwind,TN,TL]=getParam(Qp);
TUV=getData(Qd);

%Dt=30000;
h=TUV(2,1)-TUV(1,1);
% h=0.04;
Nt=round(Dt/h);
TUV=TUV(1:Nt,:);

TT=TUV(:,1); % time 
U=TUV(:,2); % original signal
V=TUV(:,3); % predicted signal

% Computing maxima of the original signal
TFU=islocalmax(U);
TUm=TT(TFU);
Um=U(TFU);
TUm1=TUm(Um>0);

% Computing maxima of the  predicted signal
TFV=islocalmax(V);
TVm=TT(TFV);
Vm=V(TFV);
TVm1=TVm(Vm>0);

figure (20)
X1a=[37 0.9];
Xat=X1a(1);
Yat=X1a(2);

X1b=[-0.95 0.9];
Xbt=X1b(1);
Ybt=X1b(2);

X1c=[-0.95 0.9];

Xct=X1c(1);
Yct=X1c(2);

X1d=[40 60];
Xdt=X1d(1);
Ydt=X1d(2);

% Short term prediction
%TN=800;
NN=round(TN/h);
%NN=20000;
subplot(3,2,[1 2]);
wideAx=gca;
axes(wideAx);
hold on
%plot(TT(1:NN),U(1:NN),'-', 'color',color_ex)
%plot(TT(1:NN),V(1:NN),'-', 'color',color_rc)

plot(TT(1:NN),U(1:NN),'-','color',color_ex,'LineWidth',lw2)
plot(TT(1:NN),V(1:NN),'-','color',color_rc,'LineWidth',lw2)

hold off
ylim([-1.1 1.1]);
xlim([0,TN]);
xlabel('time [ms]','Interpreter','latex','FontSize',10);
ylabel('$u$, $\nu$','Interpreter','latex','FontSize',10);
text(Xat,Yat,'(a)','FontSize',10);
leg1=legend('original', 'reservoir', 'Orientation', 'horizontal','Interpreter','latex','FontSize',10);
leg1.Position = [0.49 0.95 0.1 0.01];
set(wideAx,'XScale','linear');
set(wideAx,'YScale','linear');

set(wideAx,'XTick',[0 200 400 600 800]);
set(wideAx,'XTickLabel',{'0','200','400','600','800','Interpreter','latex','FontSize',10})

set(wideAx,'YTick',[-1 0 1]) 
set(wideAx,'YTickLabel',{'-1','0','1','Interpreter','latex','FontSize',10})


box on

% Long term prediction
%h=0.04;
%TL=5000;
NL=round(TL/h);
tau1=0.6;
Ntau=round(tau1/h);

subplot(3,2,3)
ha=gca;
axes(ha);
hold on
%-------------------------------------------------------------
%N=size(U(1:NL-Ntau),1);
N=NL-Ntau;
%segLength = 80;
x1=U(Ntau+1:NL);
x2=U(1:NL-Ntau);

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
subplot(3,2,4)
ha=gca;
axes(ha);
hold on;
%-------------------------------------------------------------
%N=size(U(1:NL-Ntau),1);
N=NL-Ntau;
%segLength = 80;
x1=V(Ntau+1:NL);
x2=V(1:NL-Ntau);

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
% Plot a map of interspike intervals for original and predicted signals
subplot(3,2,5)
ha=gca;
axes(ha);
plot(TUm1(2:end-1)-TUm1(1:end-2),TUm1(3:end)-TUm1(2:end-1),'.','color',color_ex,'MarkerSize',lw3)
hold on
plot(TVm1(2:end-1)-TVm1(1:end-2),TVm1(3:end)-TVm1(2:end-1),'.', 'color',color_rc,'MarkerSize',lw3)
hold off;
xlim([20, 65])
ylim([20, 65])

set(ha,'XScale','linear');
set(ha,'YScale','linear');
xlabel('$T_{n}$ [ms]','Interpreter','latex','FontSize',10);
ylabel('$T_{n+1}$ [ms]','Interpreter','latex','FontSize',10);
text(Xdt,Ydt,'(d)','FontSize',10);

set(ha,'XTick',[20  40  60]);
set(ha,'XTickLabel',{'20','40','60','Interpreter','latex','FontSize',10})
set(ha,'YTick',[20 40 60]) 
set(ha,'YTickLabel',{'20','40','60','Interpreter','latex','FontSize',10})
%-----------------------------------------------------------------------------------------------
box on
end 

