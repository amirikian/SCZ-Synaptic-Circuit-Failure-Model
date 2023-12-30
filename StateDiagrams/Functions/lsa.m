function ls = lsa(lambda, omega, mf, sp)
%   Computes quantities for LS analysis
%   
%     Input:
%         lambda - instability growth rate
%         omega  - oscillation frequency
%         mf     - structure providing MF synaptic currents and response function slope
%         sp     - structure providing synaptic time constants
%     Output:
%         ls     - structure providing quantities for LS analysis

global E I

ls.SS_AMPA = exp(-lambda*sp.tauD(E)*1e-3)/sqrt(((1+lambda*sp.tauAr*1e-3)^2+(omega*sp.tauAr*1e-3)^2)*((1+lambda*sp.tauAd*1e-3)^2+(omega*sp.tauAd*1e-3)^2));
ls.SS_NMDA = exp(-lambda*sp.tauD(E)*1e-3)/sqrt(((1+lambda*sp.tauNr*1e-3)^2+(omega*sp.tauNr*1e-3)^2)*((1+lambda*sp.tauNd*1e-3)^2+(omega*sp.tauNd*1e-3)^2));
ls.SS_GABA = exp(-lambda*sp.tauD(I)*1e-3)/sqrt(((1+lambda*sp.tauGr*1e-3)^2+(omega*sp.tauGr*1e-3)^2)*((1+lambda*sp.tauGd*1e-3)^2+(omega*sp.tauGd*1e-3)^2));

ls.phiAMPA = omega*sp.tauD(E)*1e-3 + atan2(omega*sp.tauAr*1e-3,(1+lambda*sp.tauAr*1e-3)) + atan2(omega*sp.tauAd*1e-3,(1+lambda*sp.tauAd*1e-3));
ls.phiNMDA = omega*sp.tauD(E)*1e-3 + atan2(omega*sp.tauNr*1e-3,(1+lambda*sp.tauNr*1e-3)) + atan2(omega*sp.tauNd*1e-3,(1+lambda*sp.tauNd*1e-3));
ls.phiGABA = omega*sp.tauD(I)*1e-3 + atan2(omega*sp.tauGr*1e-3,(1+lambda*sp.tauGr*1e-3)) + atan2(omega*sp.tauGd*1e-3,(1+lambda*sp.tauGd*1e-3));

ls.xxEA = mf.AA(E)*ls.SS_AMPA*mf.iAMPA(E)/mf.iSYN(E);
ls.xxEN = mf.AA(E)*ls.SS_NMDA*mf.iNMDA(E)/mf.iSYN(E);
ls.xxIG =-mf.AA(I)*ls.SS_GABA*mf.iGABA(I)/mf.iSYN(I);
end