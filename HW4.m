function HW4
clc; clear;
p3
end

function p1
r1 = 1;
r2 = 5.2;

%% Part A
%syms th real
%c = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th));
%s = (r1 + r2 + c)/2;
%am = vpa(s/2);

%% Part B
mu = 1.327e11; % Standard Gravitational Parameter of the Sun [km^3/s^2]
%mu = 1; % Standard Gravitational Parameter of the Sun [AU^3/TU^2]
th = 150;

c = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th));
s = (r1 + r2 + c)/2;
a_m = s/2;
fprintf('a_m   : %f AU\n',a_m);

% Calculate t_m
alpha_m = 2*asin(sqrt(s/(2*a_m)));
beta_m = 2*asin((sqrt((s-c)/(2*a_m))));
t_m = ((sqrt((au2km(s))^3/8))*(pi-beta_m + sin(beta_m)))/sqrt(mu);
fprintf('t_m   : %f years\n',t_m/31536000);

% Calculate t_f
a=5; % Semi-major Axis [au]
alpha_o = 2*asin(sqrt(s/(2*a)));
beta_o = 2*asin((sqrt((s-c)/(2*a))));
t_f = (sqrt((au2km(a)^3)/mu)*(alpha_o-beta_o - sin(alpha_o) + sin(beta_o)));
fprintf('t_f   : %f years\n',t_f/(60*60*24*365.25));

% Calculate t_f^#
alpha_1 = 2*pi - alpha_o;
beta_1 = beta_o;
t_fup = (sqrt((au2km(a)^3)/mu)*(alpha_1-beta_1 - sin(alpha_1) + sin(beta_1)));
fprintf('t_f^# : %f years\n',t_fup/31536000);

% Calculate t_p 
t_p = (sqrt(2)/(3*sqrt(mu)))*(au2km(s)^(3/2) - sign(sind(th))*(au2km(s)-au2km(c))^(3/2));
fprintf('t_p   : %f years\n',t_p/31536000)

%% Part C

gamma = asind((r2*sind(th))/c);

ui = [1;0];
uc = [-cosd(gamma); sind(gamma)];

A_o = sqrt(1/(4*a))*cot(alpha_o/2);
B_o = sqrt(1/(4*a))*cot(beta_o/2);

v1_o = (B_o+A_o)*uc + (B_o-A_o)*ui;
fprintf('v_1   : %f i + %f j AU/TU \n',v1_o)

A_1 = sqrt(1/(4*a))*cot(alpha_1/2);
B_1 = sqrt(1/(4*a))*cot(beta_1/2);

v1_1 = (B_1+A_1)*uc + (B_1-A_1)*ui;
fprintf('v_1^# : %f i + %f j AU/TU \n',v1_1)

%% Part D
fprintf('|v_1|   : %f AU/TU \n',norm(v1_o))
fprintf('|v_1^#| : %f AU/TU \n',norm(v1_1))

%% Part E
p = ((4*a*(s-r1)*(s-r2))/(c^2))*(sin((alpha_o+beta_o)/2))^2;
p_tild = ((4*a*(s-r1)*(s-r2))/(c^2))*(sin((alpha_1+beta_1)/2))^2;
fprintf('p       : %f AU \n',p)
fprintf('p_tilda : %f AU \n',p_tild)

e = sqrt(1-(p/a));
e_tild = sqrt(1-(p_tild/a));
fprintf('e       : %f \n',e)
fprintf('e_tilda : %f \n',e_tild)
end

function p2
r1 = 1; %Earth [au]
a_mars = 1.5237; % Mars [au]
e_mars = 0.0934;
f = 210;
r2 = (a_mars*(1-e_mars^2))/(1+e_mars*cosd(f));
th = 120; % deg
mu = 1;

%% Part A
c = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th));
s = (r1 + r2 + c)/2;
a_m = s/2;

mu = 1.327e11; % Standard Gravitational Parameter of the Sun [km^3/s^2]
%mu = 1; % Standard Gravitational Parameter of the Sun [AU^3/TU^2]

% Calculate t_m
beta_m = 2*asin((sqrt((s-c)/(2*a_m))));
t_m = ((sqrt((au2km(s)^3)/8))*(pi-beta_m + sin(beta_m)))/sqrt(mu);


%% Part B
p = ((4*a_m*(s-r1)*(s-r2))/(c^2))*(sin((pi+beta_m)/2))^2;
e = sqrt(1-(p/a_m));


%% Part C
T_mars = 2*pi*sqrt((au2km(a_mars)^3)/mu); %Mars Orbital Period [s]

syms E2 M tx f1
f2 = 210;
E2 = vpa(solve(tand(f2/2) == sqrt((1+e_mars)/(1-e_mars))*tan(E2/2),E2));
M2 = solve(M == E2 - e_mars*(sin(E2)),M);
tx = solve(M2 == sqrt(mu/au2km(a_mars)^3)*(t_m+tx),tx); %Time of mars at departure past perigee


M1 = sqrt(mu/au2km(a_mars)^3)*tx;
[iterations, values] = KeplersEqnNewtonMethod(M1,e_mars,M1,10^-14);
E1 = values(end);
f1 = rad2deg(double(vpa(solve(tan(f1/2) == sqrt((1+e_mars)/(1-e_mars))*tan(E1/2),f1))));

b_mars = a_mars*sqrt(1-e^2);

ra = [a_mars*cos(E1)-e_mars; b_mars*sin(E1);0];

%OR

%Define a new elapsed time past epoch.
t_apogee = T_mars/2; %Time at apogee
t_o = t_apogee-(T_mars + tx);

a = au2km(a_mars);
r_o = [a*(1+e_mars); 0; 0];
v_o = [0; sqrt(mu*((2/norm(r_o))-(1/a))); 0];
H_o = mu/norm(r_o)^3;
P_o = 0;
rb = km2au(r_o*(1-(t_o^2/2)*H_o) + v_o*(t_o - (t_o^3/6)*H_o));



%% Part D

gamma = asind((r2*sind(th))/c);
ui = [0;1];
uc = [-cosd(gamma); -sind(gamma)]
A = sqrt(1/(4*a_m))*cot(pi/2);
B = sqrt(1/(4*a_m))*cot(beta_m/2);
v1 = (B+A)*uc + (B-A)*ui;
%% Part E
deltaV1 = v1 - [-1;0];

%% DISPLAYZ 
fprintf('t_m    : %f days\n',t_m/(60*60*24));
fprintf('e      : %f \n',e)
fprintf('T_mars : %f days\n\n',T_mars/(60*60*24));
disp('Mars Departure Position via Potision Equation')
fprintf('tx       : %f days\n',tx/(60*60*24));
fprintf('r_mars,d : %f i + %f j + %f z AU \n\n',ra)
disp('Mars Departure Position Via F&G Functions')
fprintf('t0       : %f days\n',t_o/(60*60*24));
fprintf('r_mars,d : %f i + %f j + %f z AU \n\n',rb)
disp('Part D')
fprintf('v_1   : %f i + %f j AU/TU \n',v1)
fprintf('|v_1| : %f AU/TU\n\n',norm(v1))
disp('Part E')
fprintf('Delta V_1   : %f AU/TU \n',norm(deltaV1))
fprintf('Delta V_1   : %f km/s \n',norm(deltaV1)*29.8)
end

function p3
r1 = 9.538; % Saturn [au]
r2 = 1.523; % Mars [au]
th = 90; % deg

%% Part A
c = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th));
s = (r1 + r2 + c)/2;
a_m = s/2;
fprintf('a_m   : %f AU\n',a_m);
% mu = 1.327e11; % Standard Gravitational Parameter of the Sun [km^3/s^2]
 mu = 1; % Standard Gravitational Parameter of the Sun [AU^3/TU^2]

 %% Part B
% Calculate t_m
alpha_m = 2*asin(sqrt(s/(2*a_m)));
beta_m = 2*asin((sqrt((s-c)/(2*a_m))));
t_m = ((sqrt((s)^3/8))*(pi-beta_m + sin(beta_m)))/sqrt(mu);
fprintf('t_m   : %f TU\n',t_m);

%% Part C
% Calculate t_p 
t_p = (sqrt(2)/(3*sqrt(mu)))*(s^(3/2) - sign(sind(th))*(s-c)^(3/2));
fprintf('t_p   : %f TU\n',t_p)

%% Part D
tf = 67.12; % TU
a=fzero(@lambert,10);
fprintf('a       : %f AU\n',a)

%% Part F
alpha = 2*pi - real(2*asin(sqrt(s/(2*a))));
beta = 2*asin((sqrt((s-c)/(2*a))));
p = ((4*a*(s-r1)*(s-r2))/(c^2))*(sin((alpha+beta)/2))^2;
e = sqrt(1-(p/a));
fprintf('e       : %f \n',e)

%% Part G

gamma = asind((r2*sind(th))/c);

ui = [1;0];
uc = [-cosd(gamma); sind(gamma)];
A = sqrt(1/(4*a))*cot(alpha/2);
B = sqrt(1/(4*a))*cot(beta/2);

v1 = (B+A)*uc + (B-A)*ui;
fprintf('v_1   : %f i + %f j AU/TU \n',v1)
fprintf('v_1   : %f i + %f j km/s \n',v1*29.8)
fprintf('|v_1|   : %f AU/TU \n',norm(v1))

aH = (r1+r2)/2;
sqrt(mu*((2/r1)-(1/aH)));
sqrt(mu/r1);

deltaV1 = sqrt(mu*((2/r1)-(1/aH))) - sqrt(mu/r1);
fprintf('Delta v_1   : %f AU/TU \n',deltaV1)
fprintf('Delta v_1   : %f km/s \n',deltaV1*29.8)

%}
end


function p4
r1 = 1; %Earth [au]
r2 = 1.524; % Mars [au]
a = 1.36; % [au]
a_m = 1.14; % [au]
th1 = 107; % deg
th2 = 253; % deg
mu = 1;
%% Part A

c1 = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th1));
s1 = (r1 + r2 + c1)/2;

c2 = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th2));
s2 = (r1 + r2 + c2)/2;

c = c1;
s = s1;
mu = 1.327e11; % Standard Gravitational Parameter of the Sun [km^3/s^2]
%mu = 1; % Standard Gravitational Parameter of the Sun [AU^3/TU^2]

% Calculate t_m
beta_m = 2*asin((sqrt((s-c)/(2*a_m))));
t_m = ((sqrt((au2km(s))^3/8))*(pi-beta_m + sin(beta_m)))/sqrt(mu);
fprintf('t_m   : %f days\n',t_m/(60*60*24));

% Calculate t_f
alpha_o = 2*asin(sqrt(s/(2*a)))
beta_o = 2*asin((sqrt((s-c)/(2*a))))
t_f = ((sqrt((au2km(s))^3/8))*(alpha_o-beta_o - sin(alpha_o) + sin(beta_o)))/sqrt(mu);
fprintf('t_f   : %f days\n',t_f/(60*60*24));

% Calculate t_f^#
alpha_1 = 2*pi - alpha_o;
beta_1 = beta_o;
t_fup = ((sqrt((au2km(s))^3/8))*(alpha_1-beta_1 - sin(alpha_1) + sin(beta_1)))/sqrt(mu);
fprintf('t_f^# : %f days\n',t_fup/(60*60*24));

%}
end

function x = au2km(s)
x = s*149598000;
end

function x = km2au(s)
x = s/149598000;
end

function f = lambert(a)
r1 = 9.538; % Saturn [au]
r2 = 1.523; % Mars [au]
th = 90; % deg
c = 9.658828759223345;
s = 10.359914379611673;
alpha = 2*pi - real(2*asin(sqrt(s/(2*a))));
beta = 2*asin((sqrt((s-c)/(2*a))));

tf = 10.69;
f = tf - (a^(3/2)/(2*pi))*(alpha-beta-sin(alpha)+sin(beta));
end