%% Newton Approach
function [counter, values] = KeplersEqnNewtonMethod(M,e,E_0,tol)
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

%% Background Functions
function f = fn(M,En,e)
% Evaluate the given function.
 f = M-(En-e*sin(En));
end

function fp = fnp(M,En,e)
% Evaluate the derivative of the given function.
 fp = -1 + e*cos(En);
end