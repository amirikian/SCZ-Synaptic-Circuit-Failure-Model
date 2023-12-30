function F = root_lambda_omega(x)
%    LS equations for solving instability growth rate and oscillation frequency 
%   
%     Input: 2 unknowns
%         x(1) - instability growth rate
%         x(2) - oscillation frequency
%
%     Output: l.h.s. of 2 equations
%         F(1) - LS eq 1
%         F(2) - LS eq 2

global sp mf

lambda = x(1)*1e3;
omega = abs(x(2)*1e3);

ls = lsa(lambda, omega, mf, sp);

F(1) = ls.xxEA*cos(ls.phiAMPA) + ls.xxEN*cos(ls.phiNMDA) - ls.xxIG*cos(ls.phiGABA) - 1;
F(2) = ls.xxEA*sin(ls.phiAMPA) + ls.xxEN*sin(ls.phiNMDA) - ls.xxIG*sin(ls.phiGABA);
end