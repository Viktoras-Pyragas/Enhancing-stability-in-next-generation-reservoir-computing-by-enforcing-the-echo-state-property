function [Xnow,Xdel,Xnow_ex,Xdel_ex] = plot_MG_dyn(Xtot)

h=0.01;
%--------------------------------------------------
%Tpr1=Xtot(:,1);
U1=Xtot(:,2); % time series of original system
XP1=Xtot(:,3); % time series of autonomous NGRC system

dx=abs(XP1-U1);
indx2=find(dx>0.02,1,'first');

% Valid prediction time
Tvp=indx2*h;
disp('Valid prediction time:');
fprintf('Tvp = %.2f\n',Tvp);

tau1=2.0; % delay time of system
dN=round(tau1/h);
Dtau=400.0; % time for plot in delayed coordinates
Ndt=round(Dtau/h);

% For NGRC system:
Xnow=XP1(dN+1:Ndt,1);
Xdel=XP1(1:Ndt-dN,1);

% For original system:
Xnow_ex=U1(dN+1:Ndt,1);
Xdel_ex=U1(1:Ndt-dN,1);

end