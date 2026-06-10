function plot_bif_tot_v

load('result_ngrc_v6t.mat','XPMaxArray'); % without restriction;
%load('result_ngrc_v6.mat','XPMaxArray');
%load('result_ngrc_witout_Lc_v6.mat','XPMaxArray');
%load('result_ngrc_v6_Ndel5_deg7.mat','XPMaxArray');
%load('result_ngrc_v6_Ndel6_deg7.mat','XPMaxArray');
%load('result_ngrc_v6_Ndel7_deg7.mat','XPMaxArray');
%load('result_ngrc_v6_Ndel5_deg8_Nb51.mat','XPMaxArray');
%load('result_ngrc_v6_Ndel5_deg8_Nb101.mat','XPMaxArray');


h=0.05; % step of integration
dn=1;
bv=XPMaxArray(1:dn:end,1); % array of parameters;
Ymx=h*XPMaxArray(1:dn:end,3); % array of time intervals;

nP=5;
%pbv=linspace(0.585,0.745,nP); % with nP=7; (successfull)
pbv=linspace(3.0,3.5,nP); % with nP=7; (successfull)
Yv=[0; 100];

load('result_exact_v.mat','XPMaxArray');

bvt=XPMaxArray(:,1); % array of parameters;
Ymxt=XPMaxArray(:,3); % array of time intervals;

figure
subplot(2,1,1)
hold on
plot(bvt,Ymxt,'b*','MarkerSize',1); % exat HMR
for np=1:nP
plot([pbv(1,np); pbv(1,np)],Yv(:,1),'r--','LineWidth',1);
end
%xlim([0.2 1]);
%xlim([3.0 3.5]);
ylim([0 100]);
hold off
subplot(2,1,2)
hold on
plot(bv,Ymx,'m*','MarkerSize',1); % NGRC from HMR
for np=1:nP
plot([pbv(1,np); pbv(1,np)],Yv(:,1),'r--','LineWidth',1);
end
%xlim([0.2 1.0]);
%xlim([3.0 3.5]);
%ylim([0 100]);
hold off

end