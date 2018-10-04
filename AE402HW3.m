function AE402HW3()
format short
clc; clear;
P3()
end

function P1()
rr = [6045; 3490; 0];
vv = [-2.457; 6.618; 2.533];
mu = 3.986004418e5; %Standard Graviational Parameter km^3 s^-2

%% Calculate Semi Major Axis
r = norm(rr)
v = norm(vv)

a = r/(2-((r*v^2)/mu))

%% Calculate Eccentricity
ee = (1/mu)*(((v^2-(mu/r))*rr)-(dot(rr,vv)*vv))
e = norm(ee)
%% Calculate Incline

hh = cross(rr,vv)
h = norm(hh)

i = acos(hh(3)/h)

%% Ascending Node

nn = cross([0;0;1],hh/h)
n = norm(nn)

if nn(2) >= 0 
    omega = acos(nn(1)/n)
elseif nn(2) < 0 
    omega = 2*pi - acos(nn(1)/n)
end


%% Argument of Periapse
if ee(3) > 0
    w = acos((dot(nn,ee))/(n*e))
elseif ee(3) < 0
    w = 2*pi - acos((dot(nn,ee))/(n*e))
else
    w = 0
end

%% True Anomaly
dot(rr,vv)
if dot(rr,vv) >= 0
    f0 = acos(dot(ee,rr)/(e*r))
elseif dot(rr,vv) < 0
    f0 = 2*pi - acos(dot(ee,rr)/(e*r))
end
    
end

function P2
a = 1.1*6378;
e = 0.05;
i = 45;
RAAN = 0;
w = 20;
f = 10;
rd = 1;
mu = 3.986004418e5;

[rr,vv] = OMtoRV(mu,a,e,i,RAAN,w,f, rd)


r = norm(rr)
l  = rr(1)/r 
m = rr(2)/r 
n =rr(3)/r        % Direction cosines

delta = asind(n);                    % Declination
% Right ascension:
if (m >0)
    alpha = acosd(l/cosd(delta));
else
    alpha = 360 - acosd(l/cosd(delta));
end
fprintf('Right ascension = %4.2f [deg] \n',alpha);
fprintf('Declination = %4.2f [deg] \n',delta);
fprintf('Right ascension = %4.2f [rad] \n',deg2rad(alpha));
fprintf('Declination = %4.2f [rad] \n',deg2rad(delta));
end

function P3()
syms x1 x2 x3 real


a = 384400; % semi major axis of moon orbit in km
%{
%% Earth Craft Moon
m1 = 5.974e24; %Mass of Earth [kg]
m2 = 0; %Mass of Moon [kg]
m3 = 7.348e22;
M = m1+m2+m3;

X = solvequintic(m1,m2,m3);
syms x1 x2 x3 real
eqn1 = a == x3+x1;
eqn2 = m1*x1 + m2*x2 + m3*x3 == 0;
eqn3 = X == (x3-x2)/(x2-x1);

sol = solve([eqn1, eqn2, eqn3],[x1,x2,x3]);

x1 = sol.x1
x2 = sol.x2
x3 = sol.x3

L1 = x2-x1
%}
%% Earth Moon Craft
%{
m1 = 5.974e24; %Mass of Earth [kg]
m2 = 7.348e22; %Mass of Moon [kg]
m3 =0;
M = m1+m2+m3;;
X = solvequintic(m1,m2,m3);
x1 = -1/M*(m2*a + m3*a*(1+X))
x2 = a+x1
x3 = 1/M*(m2*a*X + m1*a*(1+X))
x = [x1;x2;x3];
L2 = x3-x2
%}
%%{
%% Craft Earth Moon
m1 = 0;
m2 = 5.974e24;
m3 = 7.348e22;
M = m1+m2+m3;;
X = solvequintic(m1,m2,m3);
syms x1 x2 x3 real
eqn1 = a == x3-x2;
eqn2 = m1*x1 + m2*x2 + m3*x3 == 0;
eqn3 = X == (x3-x2)/(x2-x1);

sol = solve([eqn1, eqn2, eqn3],[x1,x2,x3])

x1 = sol.x1
x2 = sol.x2
x3 = sol.x3

L3 = x1+x2

%}
end

function P4()
syms x1 x2 x3 real
m1 = 1.989e30; %Mass of sun [kg]
m2 = 5.974e24; %Mass of earth [kg]
m3 =0;
a = 1; % [Au]
M = m1+m2+m3;
X = solvequintic(m1,m2,m3);


x1 = -1/M*(m2*a + m3*a*(1+X));
x3 = 1/M*(m2*a*X + m1*a*(1+X));
x2 = a+x1

L2 = x3-x2;

end

function X=solvequintic(m1,m2,m3)
    syms x
    XX=vpa(solve((m1+m2)*x.^5 + (3*m1+2*m2)*x.^4+(3*m1+m2)*x.^3 - (m2+3*m3)*x.^2 - (2*m2 + 3*m3)*x - (m2+m3)));
    XX(imag(XX) ~= 0) = [];
    
    X=XX;
end