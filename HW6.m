function HW6
clc; clear;
P2
end

function P1
r1 = au2km(1);
r2 = au2km(5.2);
mu_sun = 1.327e11;
mu_earth = 3.986e5;
mu_jup = 126686000;

[a_H,t_H,dv_1, dv_2, dv_tot] = HohmannTransfer(r1,r2,mu_sun);

v_earth = sqrt(mu_sun/r1)
v_jup = sqrt(mu_sun/r2)
Vinfe = dv_1;
vinfj = dv_2;

rp=74491;
vp = sqrt(vinfj^2 + ((2*mu_jup)/rp))

delta = 2*asind(1/(1+((rp*vinfj^2)/mu_jup)))

vsco = sqrt(vinfj^2 + v_jup^2 - 2*vinfj*v_jup*cosd(delta))

gamma = asind((vinfj/vsco)*sind(delta))

f = 2*gamma

%h = r2*vsco*cosd(gamma)
h = km2au(r2)*vsco/29.8 * cosd(gamma)

p = h^2
r_sat = 9.5388;
f2 = acosd((p/r_sat)-1)

f2 - f
end

function P2
%% Part A
T = (2/3)*365*24*60*60;
mu_sun = 1.327e11;
syms a e real
a = vpa(solve(T == 2*pi*sqrt((a^3)/mu_sun)));

fprintf('a : %f km \na : %f AU \n',a, km2au(a));
rp = au2km(1)-a;
e = 1-(rp/a);
fprintf('e : %f \n',e);

%% Part B
mu_earth = 3.986e5;

V_earth = 29.78;

V_ah = sqrt(mu_earth*((2/(a*(1-e)))-(1/a)))

end

function P3
r1 = 1;
r2 = 0.723;
mu = 1;

v1 = sqrt(mu*((2/r1)-(1/1.1)))
v1 = sqrt(mu*((2/r2)-(1/1.1)))
end

function x = au2km(s)
x = s*149598000;
end

function x = km2au(s)
x = s/149598000;
end