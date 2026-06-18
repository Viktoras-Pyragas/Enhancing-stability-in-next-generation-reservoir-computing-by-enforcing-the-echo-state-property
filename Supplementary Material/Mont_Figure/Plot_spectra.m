function Plot_spectra(Qp,Qd)
% - Computing power spectra of the original and recunstructed signals
[lw2,lw3,lwd,color_ex,color_rc,a1,a2,segLength,Dt,Nwind]=getParam(Qp);
TUV=getData(Qd);
h=TUV(2,1)-TUV(1,1);
Nt=round(Dt/h);
TUV=TUV(1:Nt,:);

U=TUV(:,2); % original signal
V=TUV(:,3); % predicted signal
t=TUV(:,1); % time

fs=1/(t(2)-t(1));

X1e=[0.3 10];
Xet=X1e(1);
Yet=X1e(2);

figure (20)

subplot(3,2,6)
ha=gca;
[PowerU, FrequencyU] = WelchPowerSpectralDensity(U, [], Hann(round(length(t) / Nwind)), 0., fs); 
semilogy(gca, FrequencyU, (abs(PowerU)), '.', 'color',color_ex, 'LineWidth', lw2)
hold on;
   [PowerV, FrequencyV] = WelchPowerSpectralDensity(V, [], Hann(round(length(t) / Nwind)), 0., fs); 
   semilogy(gca, FrequencyV, (abs(PowerV)), '.', 'color',color_rc,'LineWidth', lw2)
   hold off;
  xlabel('Frequency [kHz]','Interpreter','latex','FontSize',10);
  ylabel('Power','Interpreter','latex','FontSize',10);
xlim([0 2])
ylim([1e-9,1e+2])

set(ha,'XTick',[0  1  2]);
set(ha,'XTickLabel',{'0','1','2','Interpreter','latex','FontSize',10})
set(ha,'XScale','linear');

set(ha,'YTick',[1e-8 1e-5 1e-2 10]) 
set(ha,'YTickLabel',{'10^{-8}','10^{-5}','10^{-2}','10^{1}','Interpreter','latex','FontSize',10})
set(ha,'YScale','log');

text(Xet,Yet,'(e)','FontSize',10);    
%-----------------------------------------------------------------------------------------------
box on    
end % function ends