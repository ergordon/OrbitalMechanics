format long g

rr = [6045; 3490; 0];
vv = [-2.457; 6.618; 2.533];
mu = 3.986004418e5; %Standard Graviational Parameter m^3 s^-2

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
    
RVtoOM(rr,vv, mu, 0);