function [a_H,t_H,dv_1, dv_2, dv_tot] = HohmannTransfer(r1,r2,mu)
a_H = (r1+r2)/2;
t_H = pi*sqrt((a_H^3)/mu);

dv_1 = sqrt(mu*((2/r1)-(1/a_H))) - sqrt(mu/r1);
dv_2 = sqrt(mu/r2) - sqrt(mu*((2/r2)-(1/a_H)));
dv_tot = abs(dv_1) + abs(dv_2);

disp('Hohmann Transfer')
fprintf('a_H : %f km\n',a_H);
fprintf('t_H : %f days\n',t_H/(60*60*24));
fprintf('dV1 : %f km/s\n',dv_1);
fprintf('dV2 : %f km/s\n',dv_2);
fprintf('dVtot : %f km/s\n\n\n',dv_tot);
end