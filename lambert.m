function [f,c,s,am,alpha,beta] = lambert(a, r1, r2, th,tf)
c = sqrt(r1^2 +r2^2 - 2*r1*r2*cosd(th))
s = (r1 + r2 + c)/2
am = s/2
alpha = 2*asin(sqrt(s/(2*a)))
beta = 2*asin((sqrt((s-c)/(2*a))))

f = tf - ((a^1.5)/(2*pi))*(alpha-beta-sin(alpha)+sin(beta));
end