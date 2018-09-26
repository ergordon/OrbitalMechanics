clc; clear; 
format long g
%{
An object is detected by AF Space Command in Colorado Springs. Its position
and velocity are determined below to be:
%}

% Given
r = [-4743;4743;0];
v = [-5.879;-4.223;0];
mu_E = 3.9860044188e5;

% a) What is the semi-major axis of the object's orbit?
a = 1/((-((norm(v))^2/mu_E))+(2/norm(r)))

% b) What is the eccentricity of the orbit?
e = sqrt(1-((norm(cross(r,v)))^2/(mu_E*a)))

% c.a) What is the instantaneous true anomaly, f?
f = acos((((a*(1-e^2))/norm(r))-1)/e)
r_dot = dot(r,v)/norm(r) %Positive therefore f is positive

% c.b) What is the Instantaneous Mean Anomaly, M?
E = acos((1-(norm(r)/a))/e)
M = E - e*sin(E)

% d) How much time (in minutes) will pass between detection of the object
% and its impact on the earths surface)
r_impact = 6378
n = sqrt(mu_E/a^3)

f_impact = acos((((a*(1-e^2))/r_impact)-1)/e)
E_impact = 2*atan(sqrt((1-e)/(1+e))*tan(f_impact/2))
M_impact = E_impact - e*sin(E_impact)

t = M_impact/n
t_min = (t/60)

% e) What will the speed of the object be at impact?

v_impact = sqrt(mu_E*((2/r_impact)-(1/a)))