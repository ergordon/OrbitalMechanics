%% Laguerre Approach
function [counter, values] = KeplersEqnLaguerreMethod(M,e,E_0,tol)
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
 fpp = -e*sin(En);
end