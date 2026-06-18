load('Xtot_record.mat');

tau=0.6;
h=0.01;
Ntau=round(tau/h);

Xd=Xtot(1+Ntau:end,2);
Xn=Xtot(1:end-Ntau,2);

Dtau=150;
Nd=round(Dtau/h);

Xn=Xn(1:Nd,1);
Xd=Xd(1:Nd,1);

figure
hold on
plot(Xn,Xd,'b-');
hold off