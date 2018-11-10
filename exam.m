function exam
clc; clear
gamma = 1.25
Tc = 2500
R = 8314
MW = 20.5
P = 7.07e6

cstar = sqrt((1/gamma)*((2/(gamma+1))^((gamma+1)/(gamma-1)))* ((R*Tc)/(MW)))

syms A PP real
n = .7
a = 3.752e-6
rho = 1900

eqn1 = solve(PP == vpa((A*((a*rho)/(sqrt((gamma/(R*Tc))*(2/(gamma+1))^((gamma+1)/(gamma-1))))))^(1/(1-n))),A)

integrate(@(PP) 0.00607/(PP^(0.7)),P,2*P)

Arat = vpa(solve(P == (A*((a*rho)/(sqrt((gamma/(R*Tc))*(2/(gamma+1))^((gamma+1)/(gamma-1))))))^(1/(1-n)),A))
end
function P = Pfunc

end
