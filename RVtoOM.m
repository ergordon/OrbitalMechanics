function [energy, E, f, i, Omega, w, p, a] = RVtoOM(r,v, mu, rd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert R and V vectors into the classical orbital mechanical elements.
%
% This function takes in the following inputs:
%         r = position vector in orbit
%         v = velocity vector of object
%         mu = Standard Graviational Parameter [km^3/s^2]
%         rd = [boolean]
%                 0 = Radians
%                 1 = Degrees
% 
% This function takes in the following:
%    energy = The energy of the orbit
%         E = Eccentricity
%         f = True Anomaly at Epoch T0
%         i = Inclination
%     Omega = Longitude of Ascending Node
%         w = Argument of Periapse
%         p = Distance of Perigee
%         a = Semi Major Axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rd <= 0
    %finds h via cross product, then finds n
    h = cross(r,v);
    n = cross([0,0,1],h);

    %finds e vector and scalar. 
    e = ((norm(v)^2 - (mu/norm(r))).*r - (dot(r,v).*v))./mu;
    E = norm(e);

    %visa viva
    energy = (((norm(v).^2)./(2)) - (mu./norm(r)));

    if E ~= 1;
        a = -mu./(2*energy);
        p = a.*(1-E.^2);
    else
        p = norm(h).^2./mu;
        a = inf;
    end 

    %finds i, omega_big, omega_small, and f
    i = acos(h(3)./norm(h));
    Omega = acos(n(1)./norm(n));
    w = acos(dot(n,e)./(norm(n).*norm(e)));
    f = acos(dot(e,r)./(norm(e).*norm(r)));

    %finds conditional values 
    if n(2) < 0 ;
        Omega = 2*pi - Omega;
    end 

    if e(3) < 0 ;
        w = 2*pi - w;
    end 

    if dot(r,v) < 0 ;
        f = 2*pi - f;
    end 

    fprintf('Classical elements\nenergy = %f\n',energy)
    fprintf('ecentricity = %f\n',E)
    fprintf('true anamoly = %f\n',f)
    fprintf('i = %f\n',i)
    fprintf('Omega = %f\n',Omega)
    fprintf('w = %f\n',w)
    fprintf('p = %f\n',p)
    fprintf('a = %f\n',a)
    fprintf('In Radians')
    
else
    %finds h via cross product, then finds n
    h = cross(r,v);
    n = cross([0,0,1],h);

    %finds e vector and scalar. 
    e = ((norm(v)^2 - (mu/norm(r))).*r - (dot(r,v).*v))./mu;
    E = norm(e);

    %visa viva
    energy = (((norm(v).^2)./(2)) - (mu./norm(r)));

    if E ~= 1;
        a = -mu./(2*energy);
        p = a.*(1-E.^2);
    else
        p = norm(h).^2./mu;
        a = inf ;
    end 

    %finds i, omega_big, omega_small, and f
    i = acosd(h(3)./norm(h));
    Omega = acosd(n(1)./norm(n));
    w = acosd(dot(n,e)./(norm(n).*norm(e)));
    f = acosd(dot(e,r)./(norm(e).*norm(r)));

    %finds conditional values 
    if n(2) < 0 ;
        Omega = 360 - Omega;
    end 

    if e(3) < 0 ;
        w = 360 - w;
    end 

    if dot(r,v) < 0 ;
        f = 360 - f;
    end 

    fprintf('Classical elements\nenergy = %f\n',energy)
    fprintf('ecentricity = %f\n',E)
    fprintf('true anamoly = %f\n',f)
    fprintf('i = %f\n',i)
    fprintf('Omega = %f\n',Omega)
    fprintf('w = %f\n',w)
    fprintf('p = %f\n',p)
    fprintf('a = %f\n',a)
    fprintf('In Degrees')
end
end
