function F = root_gALL_omega( x )
%    MF & LS equations with constraints for solving synaptic conductances  
%    and oscillation frequency
%   
%    Case SYNC = 0, MF solution for 8 synaptic conductances
%     Input: 8 unknown conductances g
%         x(1) - g(rEA) recurrent AMPAR on E neuron
%         x(2) - g(rEG) recurrent GABAR on E neuron
%         x(3) - g(rEN) recurrent NMDAR on E neuron
%         x(4) - g(xEA) external  AMPAR on E neuron
%         x(5) - g(rIA) recurrent AMPAR on I neuron
%         x(6) - g(rIG) recurrent GABAR on I neuron
%         x(7) - g(rIN) recurrent NMDAR on I neuron
%         x(8) - g(xIA) external  AMPAR on I neuron
%     Output: l.h.s. of 8 equations
%         F(1) - MF eq 1
%         F(2) - MF eq 2
%         F(3) - constraint 1, equal NMDA/GABA current balance for E and I neurons
%         F(4) - constraint 2, fixed NMDA/GABA current balance
%         F(5) - constraint 3, equal AMPA/GABA current balance for E and I neurons
%         F(6) - constraint 4, fixed AMPA/GABA current balance
%         F(7) - constraint 5, equal xAMPA/GABA current balance for E and I neurons
%         F(8) - constraint 6, fixed xAMPA/'firing threshold' current balance 
%
%    Case SYNC = 1, MF & LS solution for 8 conductances and frequency at critical point lambda = 0
%     Input: 9 unknowns: 8 conductances g and oscillation frequency omega
%         x(1) - g(rEA) recurrent AMPAR on E neuron
%         x(2) - g(rEG) recurrent GABAR on E neuron
%         x(3) - g(rEN) recurrent NMDAR on E neuron
%         x(4) - g(xEA) external  AMPAR on E neuron
%         x(5) - g(rIA) recurrent AMPAR on I neuron
%         x(6) - g(rIG) recurrent GABAR on I neuron
%         x(7) - g(rIN) recurrent NMDAR on I neuron
%         x(8) - g(xIA) external  AMPAR on I neuron
%         x(9) - omega oscillation frequency
%     Output: l.h.s. of 9 equations
%         F(1) - MF eq 1
%         F(2) - MF eq 2
%         F(3) - constraint 1, equal NMDA/GABA current balance for E and I neurons
%         F(4) - constraint 2, fixed NMDA/GABA current balance
%         F(5) - constraint 3, equal AMPA/GABA current balance for E and I neurons
%         F(6) - constraint 4, fixed AMPA/GABA current balance
%         F(7) - constraint 5, equal xAMPA/GABA current balance for E and I neurons
%         F(8) - LS eq 1 
%         F(9) - LS eq 2 

global sp E I
global iAiG iNiG
global iThrE iXiT SYNC

x = abs(x);
if SYNC
    sp.g = x(1:end-1);
    omega = x(end)*1e3;
    lambda = 0;
else
    sp.g = x;
end

[nu0, mf] = mfa(sp);

if SYNC
    ls = lsa(lambda, omega, mf, sp);
end

F(1) = nu0(E) - sp.nu(E);
F(2) = nu0(I) - sp.nu(I);

F(3) = mf.iNMDA(E)/mf.iGABA(E) - mf.iNMDA(I)/mf.iGABA(I);
F(4) =-mf.iNMDA(E)/mf.iGABA(E) - iNiG;

F(5) = mf.iAMPA(E)/mf.iGABA(E) - mf.iAMPA(I)/mf.iGABA(I);
F(6) =-mf.iAMPA(E)/mf.iGABA(E) - iAiG;

F(7) = mf.ixAMPA(E)/mf.iGABA(E) - mf.ixAMPA(I)/mf.iGABA(I);

if SYNC
    F(8) = ls.xxEA*cos(ls.phiAMPA) + ls.xxEN*cos(ls.phiNMDA) - ls.xxIG*cos(ls.phiGABA) - 1;
    F(9) = ls.xxEA*sin(ls.phiAMPA) + ls.xxEN*sin(ls.phiNMDA) - ls.xxIG*sin(ls.phiGABA);
else
    F(8) = mf.ixAMPA(E)/iThrE - iXiT;
end
end