%{
A satellite is in an elliptic orbit with e=0.5. Using a starting value of
E_0 = M, determine the value of E to an accuracy of 10^-4 radians at the
time the satellite is one-quarter period past periapse passage. List all
iterations including the value of F(E) by:
%}

function HW2P1
clc; clear;
mu = 3.9860044188e5;
e = .5;
M = pi/2;
%e = 0.999;
%M = pi/7;
tol = 10^-4;

E_0 = M;
% a) Newton Method
[iterationsNewtona, valuesa] = newton(M,e,E_0,tol);
HW2P1a = [iterationsNewtona, valuesa]

% b) Laguerre Method
[iterationsLaguerreb, valuesb] = laguerre(M,e,E_0,tol);
HW2P1b = [iterationsLaguerreb, valuesb]

E_0 = eqn216(M,e)
% c.a) Newton Method
[iterationsNewtonc, valuesca] = newton(M,e,E_0,tol);
HW2P1ca = [iterationsNewtonc, valuesca]

% c.b) Laguerre Method
[iterationsLaguerrec, valuescb] = laguerre(M,e,E_0,tol);
HW2P1cb = [iterationsLaguerrec, valuescb]

% d) Solve for the True Anomaly, f.
E = valuescb(end)
f = 2*atan2(tan(E/2),sqrt((1-e)/(1+e)))
end

%% Newton Approach
function [counter, values] = newton(M,e,E_0,tol)
x = E_0;
i=0;
counter = [i];
values = [x];
    while abs(fn(M,x,e))>tol
        dx = -fn(M,x,e)/fnp(M,x,e);
        i=i+1;
        x = x+dx;
        
        counter = [counter; i];
        values = [values;x];
    end
end

%% Laguerre Approach
function [counter, values] = laguerre(M,e,E_0,tol)
x = E_0;
i=0;
n = 4;
counter = [i];
values = [x];
    while abs(fn(M,x,e))>tol
        %dx = -fn(M,x,e)/fnp(M,x,e);
        if fnp(M,x,e)>0
            dx = -(n*fn(M,x,e))/(fnp(M,x,e)+sqrt((n-1)^2*(fnp(M,x,e))^2 - n*(n-1)*fn(M,x,e)*fnpp(M,x,e)));
            i=i+1;
            x = x+dx;
        elseif fnp(M,x,e)<0
            dx = -(n*fn(M,x,e))/(fnp(M,x,e)-sqrt((n-1)^2*(fnp(M,x,e))^2 - n*(n-1)*fn(M,x,e)*fnpp(M,x,e)));
            i=i+1;
            x = x+dx; 
        else
            dx = 0;
            i=i+1;
            x = x+dx; 
        end
        counter = [counter; i];
        values = [values;x];
    end
end

%% Background Functions

function f = fn(M,En,e)
% Evaluate the given function.
 f = M-(En-e*sin(En));
end

function fp = fnp(M,En,e)
% Evaluate the derivative of the given function.
 fp = -1 + e*cos(En);
end

function fpp = fnpp(M,En,e)
% Evaluate the derivative of the given function.
 fpp = -e*cos(En);
end

function E_0 = eqn216(M,e)
% Equation 2.16 in the Orbital Mechanics Textbook
u = M+e;
E_0 = (M*(1-sin(u)+u*sin(M)))/(1+sin(M)-sin(u));
end