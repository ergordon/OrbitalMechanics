function X=solveQuintic(m1,m2,m3)
    syms x
    XX=vpa(solve((m1+m2)*x.^5 + (3*m1+2*m2)*x.^4+(3*m1+m2)*x.^3 - (m2+3*m3)*x.^2 - (2*m2 + 3*m3)*x - (m2+m3)));
    XX(imag(XX) ~= 0) = [];
    
    X=XX;
end