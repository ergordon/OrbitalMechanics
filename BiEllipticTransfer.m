function [a_1,a_2,dv_1,dv_i,dv_2,dv_tot,R,t1,t2,ttot] = BiEllipticTransfer(r1,r2,ri,mu)
a_1 = (r1+ri)/2;
a_2 = (ri+r2)/2;

dv_1 = sqrt(mu*((2/r1)-(1/a_1))) - sqrt(mu/r1);
dv_i = sqrt(mu*((2/ri)-(1/a_2))) - sqrt(mu*((2/ri)-(1/a_1)));
dv_2 = sqrt(mu/r2) - sqrt(mu*((2/r2)-(1/a_2)));
dv_tot = abs(dv_1) + abs(dv_2) + abs(dv_i);

R = r2/r1;

t1 = pi*sqrt((a_1^3)/mu);
t2 = pi*sqrt((a_2^3)/mu);
ttot = t1+t2;

disp('Bi-Elliptic Transfer')
fprintf('a_1 : %f km\n',a_1);
fprintf('a_2 : %f km\n',a_2);
fprintf('dV1 : %f km/s\n',dv_1);
fprintf('dVi : %f km/s\n',dv_i);
fprintf('dV2 : %f km/s\n',dv_2);
fprintf('dVtot : %f km/s\n',dv_tot);
fprintf('R : %f \n',R);
fprintf('t1 : %f days\n',t1/(60*60*24));
fprintf('t2 : %f days\n',t2/(60*60*24));
fprintf('t_total : %f days\n',ttot/(60*60*24));
end