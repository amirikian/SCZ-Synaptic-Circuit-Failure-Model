function F = root_nu( nu )
%    MF equations for solving population firing rates
%   
%     Input: 2 unknowns
%         nu(1) - E population rate
%         nu(2) - I population rate
%
%     Output: l.h.s. of 2 equations
%         F(1) - MF eq 1
%         F(2) - MF eq 2

global sp

sp.nu = abs(nu);
nuMF = mfa(sp);
F = nuMF - sp.nu;
end