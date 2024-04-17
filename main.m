% Script for 223 Design Project

%% Housekeeping
clear all; close all;

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

% Air (0.09)
% Mu = 2.139*10^-5; 
% Rhofl = 0.9718;
% kfl = 0.03024;
% Cpfl = 1008;

% Water (18.24)
% Mu = 1.002e-3; 
% Rhofl = 998;
% kfl = 0.598;
% Cpfl = 4182;

% R-134a (4.97)
% Mu = 2.142e-4; 
% Rhofl = 1226;
% kfl = 0.0856;
% Cpfl = 1408;

% Glycol (5.62)
% Mu = 0.0184; 
% Rhofl = 1114;
% kfl = 0.256;
% Cpfl = 2430;

% Hg (Fail)
% Mu = 1.534e-3; 
% Rhofl = 13534;
% kfl = 8.51533;
% Cpfl = 139.4;

% Isobutane (4.58)
% Mu = 1.51e-4; 
% Rhofl = 550.7;
% kfl = 0.0956;
% Cpfl = 2455;

% Ammonia (17.89)
% Mu = 1.519e-4; 
% Rhofl = 610.2;
% kfl = 0.4927;
% Cpfl = 4745;

% saturated water, saturated refrigerant, saturated ammonia, saturated ammonia, liquid metals
%% T(r) < 180 C , alpha < 10 micrometers
syms h egen
Pr = (Mu*Cpfl)/kfl;
SigmaAllow = (Kc/(1.12*sqrt(pi*alpha)))/3; 
U = 20; 
Re = (Rhofl*U*D)/Mu;
Nu = (h*D)/kfl == 2 + 0.6*(Re)^0.5*Pr^(1/3);
h = vpasolve(Nu);
eqn = (-(egen/k)*(R^2/6)) + Tinf + ((egen*R)/(3*h))+((egen*R^2)/(6*k)) == 160;
egen = vpasolve(eqn);
SigmaTheta = (2*a*E)/(1-v)*(1/5)*(egen*R^2)/(6*k)
if SigmaTheta < SigmaAllow
    disp('safe')
else
    disp('fail')
end


