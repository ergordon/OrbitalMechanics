function [Rvector,Vvector] = OMtoRV(mu,a,e,i,RAAN,w,f, rd)




%finds theta
theta = w + f;

%using vis viva, infind the scalar values for r,v
r = (a*(1-e^2))/(1+e*cosd(f))
v = sqrt(mu*((2/r)-(1/a)));
h = sqrt(mu*a*(1-e^2));

%using g,f functions it finds the vectors
Rvector = r.*[cosd(RAAN).*cosd(theta)-sind(RAAN).*sind(theta).*cosd(i), sind(RAAN).*cosd(theta)+cosd(RAAN).*sind(theta).*cosd(i), sind(theta).*sind(i)];
Vvector = -(mu./h).*[cosd(RAAN).*(sind(theta)+e.*sind(w))+sind(RAAN).*(cosd(theta)+e.*cosd(w)).*cosd(i), sind(RAAN).*(sind(theta)+e.*sind(w))-cosd(RAAN).*(cosd(theta)+e.*cosd(w)).*cosd(i), (-(cosd(theta)+e.*cosd(w)).*sind(i))];

%prints results
fprintf('Velocity: (%f, %f, %f)\n',Vvector);
fprintf('Radius: (%f, %f, %f)\n',Rvector);


end