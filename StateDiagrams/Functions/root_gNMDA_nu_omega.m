function F = root_gNMDA_nu_omega( gNnu )
%    MF & LS equations for solving population firing rates, 
%    NMDAR conductances, and oscillation frequency at critical point lambda=0
%   
%     Input: 4 unknowns
%         gNnu(1) - modulation strength of NMDAR conductances
%         gNnu(2) - E population rate
%         gNnu(3) - I population rate
%         gNnu(4) - oscillation frequency
%
%     Output: l.h.s. of 4 equations
%         F(1) - MF eq 1
%         F(2) - MF eq 2
%         F(3) - LS eq 1
%         F(4) - LS eq 2

global sp
global grEN0 grIN0 
global rEN rIN E I

lambda = 0;
gNnu = abs(gNnu);
sp.g(rEN) = gNnu(1)*grEN0;
sp.g(rIN) = gNnu(1)*grIN0;
sp.nu = gNnu(2:3);
omega = gNnu(end)*1e3;

[nuMF, mf] = eLifeMfa(sp);
ls = eLifeLsa(lambda,omega,mf,sp);

F(1) = nuMF(E)-sp.nu(E);
F(2) = nuMF(I)-sp.nu(I);
F(3) = ls.xxEA*cos(ls.phiAMPA) + ls.xxEN*cos(ls.phiNMDA) - ls.xxIG*cos(ls.phiGABA) - 1;
F(4) = ls.xxEA*sin(ls.phiAMPA) + ls.xxEN*sin(ls.phiNMDA) - ls.xxIG*sin(ls.phiGABA);
end