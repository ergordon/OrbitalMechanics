function HW2P3
clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 7200; % Semi-Major Axis [km]
e = 0.06; % Eccentricity 
mu_earth = 3.9860044188e5; % Standard Gravitational Parameter [km^3/s^2]

%Period of an Orbit
T_sec = 2*pi*sqrt(a^3/mu_earth); %[sec]
T_min = (T_sec/60); % [min]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 60; % [min]
t = t*60; % [sec]

%Calculate Mean Anomaly
M = sqrt(mu_earth/a^3)*t;

% Calculate Eccentric Anomaly using Newtons Method
E_0 = M;
tol = 10^-14;
[iterations, values] = KeplersEqnNewtonMethod(M,e,E_0,tol);
E = values(end);

%Caclulate True Anomaly
f = 2*atan2(tan(E/2),sqrt((1-e)/(1+e)));

%Verify Logics
T_min/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part C: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the Semi-Minor Axis
b = a*sqrt(1-e^2);

%Calculate Position Vector
r = [a*(cos(E)-e); (b/a)*(a*sin(E)); 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part D: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 60; % [min]
t = t*60; % [sec]

r_o = [a*(1-e); 0; 0];
v_o = [0; sqrt(mu_earth*((2/norm(r_o))-(1/a))); 0];
H_o = mu_earth/norm(r_o)^3;
P_o = 0;

(1-(t^2/2)*H_o)

(t - (t^3/6)*H_o)

r = r_o*(1-(t^2/2)*H_o) + v_o*(t - (t^3/6)*H_o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part E: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The agreement of the two methods do not show any similarities. The
% problem with this method is that we are expanding around a specific point
% that is really far from the point we are currently solving for. We will
% find that as t-t_o gets larger, the accuracy of our solution will
% decrease. A solution to this would be to include higher order terms to
% recover some accuracy at the expense of increased computation time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part F: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define a new elapsed time past epoch.
t_apogee = T_sec/2 %Time at apogee
t_o = t - t_apogee

r_o = [-a*(1+e); 0; 0]
v_o = [0; -sqrt(mu_earth*((2/norm(r_o))-(1/a))); 0]
H_o = mu_earth/norm(r_o)^3
P_o = 0


(1-(t_o^2/2)*H_o)
(t_o - (t_o^3/6)*H_o)

r = r_o*(1-(t_o^2/2)*H_o) + v_o*(t_o - (t_o^3/6)*H_o)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part G: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The estimate when expanding the taylor series around the apogee gives a
% much better result when comparing to the position from part c. This is
% likely due to the fact that the quatity t-t_o is much smaller than if we
% had set our t to be at periapsis. Apoapsis and  periapsis  are  good
% choices for expaonding around due to the fact that the velocity and
% position vecotrs are orthoganal, thus simplifying our expression.

end