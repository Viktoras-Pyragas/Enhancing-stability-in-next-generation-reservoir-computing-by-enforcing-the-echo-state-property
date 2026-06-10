function  [X,Y,Z]=HMR_signal(hT,NT,y0)
% Generuojami X HMR signalas ilgio NT su zingsniu hT
P.r=0.006;
P.nu=4;
P.kap=-1.6;
P.I=3.2;
 
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
T=linspace(0,hT*(NT-1),NT);
[~,Y1] = ode45(@(t,y) sistema(t,y,P), T,y0,opts);
format long
disp(Y1(end,:));
%format short
X=(Y1(:,1));
Y=(Y1(:,2));
Z=(Y1(:,3));
end
%---------------------------------------------------------------------
function dy=sistema(~,x,P)

   dy=   [x(2)-x(1)^3+3*x(1)^2-x(3)+P.I;
          1-5*x(1)^2-x(2);
          P.r*(P.nu*(x(1)-P.kap)-x(3))];
end 
      
  
