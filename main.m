% Script for 223 Design Project
%% Solid Constants
Rho = 2585; %kg/m3
Cp = 770; %J/kgC
k = 0.63; %W/mC
E = 47.29; %GPa
v = 0.253; 
a = 133.6*10^-7; %1/K
Kc = 0.480; %MPa/sqrt(m)
Tf = 461; %C
%egen = 1;
Tinf = 20; %C
D = 0.1; %m
R = D/2; %m
alpha = 10*10^-6; %m
%% Fluid Constants
Mu = 2.139*10^-5; 
Rhofl = 0.9718;
kfl = 0.03024;
Cpfl = 1008;
% saturated water, saturated refrigerant, saturated ammonia, saturated ammonia, liquid metals
%% T(r) < 180 C , alpha < 10 micrometers
syms h egen
Pr = (Mu*Cpfl)/kfl;
SigmaAllow = (Kc/(1.12*sqrt(pi*alpha)))/3; 
U = 20; 
Re = (Rhofl*U*D)/Mu;
Nu = (h*D)/kfl == 2 + 0.6*(Re)^0.5*Pr^(1/3);
h = vpasolve(Nu)
eqn = (-(egen/k)*(R^2/6)) + Tinf + ((egen*R)/(3*h))+((egen*R^2)/(6*k)) == 160;
egen = vpasolve(eqn)
SigmaTheta = (2*a*E)/(1-v)*(1/5)*(egen*R^2)/(6*k);
if SigmaTheta > SigmaAllow
    disp('safe')
else
    disp('fail')
end


