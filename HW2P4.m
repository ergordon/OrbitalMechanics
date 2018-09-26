function HW2P4
clc; clear;

% An earth-orbiting satellite has a period of 15.743 hr and a perigee
% radius of 2 earth radii

T = 15.743; % Period [hr]
mu_earth = 3.9860044188e5; % Standard Gravitational Parameter km^3/s^2
R_earth = 6371; % Radius of Earth [km]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part A: Determine the semimajor axis of the orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = T*60*60; % Period [sec]
r_perigee = 2*R_earth; % Radius at perigee [km]

% Semi-Major Axis
a = (mu_earth/((2*pi)/(T))^2)^(1/3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part A': Determine the position and velocity vectors at perigee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_o = [r_perigee; 0; 0];
v_o = [0; sqrt(mu_earth*((2/r_perigee)-(1/a))); 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part B: Determine the position and velocity vectors after perigee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 10;% Perigee Passage [hr]
t = t*60*60;% Perigee Passage [sec]

%Calculate Mean Anomaly
M = sqrt(mu_earth/a^3)*t

% Calculate Eccentricity
e = 1-(r_perigee/a)

% Calculate the Semi-Minor Axis
b = a*sqrt(1-e^2)

% Calculate Eccentric Anomaly using Newtons Method
E_0 = M;
tol = 10^-14

[iterations, values] = KeplersEqnNewtonMethod(M,e,E_0,tol)
E = values(end)

%Calculate f and g functions (E_o = 0 b/c at perigee)
f = 1 - (a/norm(r_o))*(1-cos(E))
g = t-sqrt(a^3/mu_earth)*(E-sin(E))


r = f*r_o + g*v_o

%This result is in the right quadrant at least. Since we are 10 hours into
%out orbit we know that we will be in quadrant 3 between 7.8715 hours and
%11.80725 hours. This ofcorurse only applies for a circular orbit but gives
%an estimate to the location and accuracy.

end