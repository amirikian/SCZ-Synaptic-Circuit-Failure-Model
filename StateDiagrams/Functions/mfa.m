function [nuMF,mf] = mfa(sp)
%   Computes MF population firing rates, synaptic currents, and
%   dimensionless slope of current-frequency response function
%   
%     Input:
%         sp     - structure providing neuron and synaptic parameters
%     Output:
%         nuMF   - population firing rates for E and I neurons
%         mf     - structure providing synaptic currents and
%                  dimensionless slope of current-frequency response function

global E I 
global Vmn0
global rEA rEG rEN xEA rIA rIG rIN xIA 

TG(E) = sp.g(rEG)*sp.cI(E)*sp.tauSG/sp.gm(E)*1e-3;
TG(I) = sp.g(rIG)*sp.cI(I)*sp.tauSG/sp.gm(I)*1e-3;

TA(E) = sp.g(rEA)*sp.cE(E)*sp.tauSA/sp.gm(E)*1e-3;
TA(I) = sp.g(rIA)*sp.cE(I)*sp.tauSA/sp.gm(I)*1e-3;

TX(E) = sp.g(xEA)*sp.cX(E)*sp.tauSA/sp.gm(E)*1e-3;
TX(I) = sp.g(xIA)*sp.cX(I)*sp.tauSA/sp.gm(I)*1e-3;

options = optimoptions(@fsolve,'Display','off');
[Vmn, ~, exitflag] = fsolve(@rootVmn, Vmn0, options);
if exitflag < 0
    error(['No solution for Vmn, Exitflag = ', num2str(exitflag)]);
end

J = 1 +  sp.gamma*exp(-sp.beta*Vmn);
rho2(E) = sp.beta*sp.g(rEN)*sp.cE(E)*sp.tauSN*(Vmn(E)-sp.VE)*(J(E)-1)/(sp.gm(E)*J(E)^2)*1e-3;
rho2(I) = sp.beta*sp.g(rIN)*sp.cE(I)*sp.tauSN*(Vmn(I)-sp.VE)*(J(I)-1)/(sp.gm(I)*J(I)^2)*1e-3;

rho1(E) = sp.g(rEN)*sp.cE(E)*sp.tauSN/(sp.gm(E)*J(E))*1e-3;
rho1(I) = sp.g(rIN)*sp.cE(I)*sp.tauSN/(sp.gm(I)*J(I))*1e-3;

S(E) = 1 + TX(E)*sp.nuX(E) + TA(E)*sp.nu(E) + (rho1(E)+rho2(E))*sp.nu(E) + TG(E)*sp.nu(I);
S(I) = 1 + TX(I)*sp.nuX(I) + TA(I)*sp.nu(E) + (rho1(I)+rho2(I))*sp.nu(E) + TG(I)*sp.nu(I);

tau = sp.tau./S;

mu(E) = (TX(E)*sp.nuX(E)+TA(E)*sp.nu(E)+rho1(E)*sp.nu(E))*(sp.VE-sp.VL)/S(E) + ...
          (rho2(E)*sp.nu(E)*(Vmn(E)-sp.VL)+TG(E)*sp.nu(I)*(sp.VI-sp.VL))/S(E);
mu(I) = (TX(I)*sp.nuX(I)+TA(I)*sp.nu(E)+rho1(I)*sp.nu(E))*(sp.VE-sp.VL)/S(I) + ...
          (rho2(I)*sp.nu(E)*(Vmn(I)-sp.VL)+TG(I)*sp.nu(I)*(sp.VI-sp.VL))/S(I);
vxAMPAsd(E) = sp.g(xEA)*abs(Vmn(E)-sp.VE)*sp.tauSA*sqrt(sp.cX(E)*sp.nuX(E)*tau(E)*1e-3)/(sp.gm(E)*sp.tau(E));
vxAMPAsd(I) = sp.g(xIA)*abs(Vmn(I)-sp.VE)*sp.tauSA*sqrt(sp.cX(I)*sp.nuX(I)*tau(I)*1e-3)/(sp.gm(I)*sp.tau(I));

vrAMPAsd(E) = sp.g(rEA)*abs(Vmn(E)-sp.VE)*sp.tauSA*sqrt(sp.cE(E)*sp.nu(E)*tau(E)*1e-3)/(sp.gm(E)*sp.tau(E));
vrAMPAsd(I) = sp.g(rIA)*abs(Vmn(I)-sp.VE)*sp.tauSA*sqrt(sp.cE(I)*sp.nu(E)*tau(I)*1e-3)/(sp.gm(I)*sp.tau(I));

vrNMDAsd(E) = sp.g(rEN)*abs(Vmn(E)-sp.VE)*sp.tauSN*sqrt(sp.cE(E)*sp.nu(E)*tau(E)*1e-3)/(sp.gm(E)*sp.tau(E))/J(E);
vrNMDAsd(I) = sp.g(rIN)*abs(Vmn(I)-sp.VE)*sp.tauSN*sqrt(sp.cE(I)*sp.nu(E)*tau(I)*1e-3)/(sp.gm(I)*sp.tau(I))/J(I);

vrGABAsd(E) = sp.g(rEG)*abs(Vmn(E)-sp.VI)*sp.tauSG*sqrt(sp.cI(E)*sp.nu(I)*tau(E)*1e-3)/(sp.gm(E)*sp.tau(E));
vrGABAsd(I) = sp.g(rIG)*abs(Vmn(I)-sp.VI)*sp.tauSG*sqrt(sp.cI(I)*sp.nu(I)*tau(I)*1e-3)/(sp.gm(I)*sp.tau(I));

sigma2 = vxAMPAsd.^2 + vrAMPAsd.^2 + vrNMDAsd.^2 + vrGABAsd.^2;

tauS = sigma2./(vxAMPAsd.^2/(sp.tauAr+sp.tauAd) + vrAMPAsd.^2/(sp.tauAr+sp.tauAd) + vrNMDAsd.^2/(sp.tauNr+sp.tauNd) + vrGABAsd.^2/(sp.tauGr+sp.tauGd));
sigma = sqrt(sigma2);
k2 = tauS./tau;

a = (sp.vThr-sp.VL-mu)./sigma.*(1+0.5*k2)+1.03*sqrt(k2)-0.5*k2;
b = (sp.vRst-sp.VL-mu)./sigma;
nuMF(E) = rspnsFunc(sp.trp(E), tau(E), a(E), b(E));
nuMF(I) = rspnsFunc(sp.trp(I), tau(I), a(I), b(I));

if nargout > 1
    mf.Vmn = Vmn;
    mf.Vhld = Vmn;
    mf.sxAMPA(E) = sp.cX(E)*sp.nuX(E)*sp.tauSA*1e-3;
    mf.sxAMPA(I) = sp.cX(I)*sp.nuX(I)*sp.tauSA*1e-3;
    mf.ixAMPA(E) = sp.g(xEA)*(mf.Vhld(E)-sp.VE)*mf.sxAMPA(E);
    mf.ixAMPA(I) = sp.g(xIA)*(mf.Vhld(I)-sp.VE)*mf.sxAMPA(I);
    
    mf.sAMPA(E) = sp.cE(E)*sp.nu(E)*sp.tauSA*1e-3;
    mf.sAMPA(I) = sp.cE(I)*sp.nu(E)*sp.tauSA*1e-3;
    mf.iAMPA(E) = sp.g(rEA)*(mf.Vhld(E)-sp.VE)*mf.sAMPA(E);
    mf.iAMPA(I) = sp.g(rIA)*(mf.Vhld(I)-sp.VE)*mf.sAMPA(I);
    
    mf.sNMDA(E) = sp.cE(E)*sp.nu(E)*sp.tauSN*1e-3;
    mf.sNMDA(I) = sp.cE(I)*sp.nu(E)*sp.tauSN*1e-3;
    mf.iNMDA(E) = sp.g(rEN)*(mf.Vhld(E)-sp.VE)*mf.sNMDA(E)/(1+sp.gamma*exp(-sp.beta*mf.Vhld(E)));
    mf.iNMDA(I) = sp.g(rIN)*(mf.Vhld(I)-sp.VE)*mf.sNMDA(I)/(1+sp.gamma*exp(-sp.beta*mf.Vhld(I)));
    
    mf.sGABA(E) = sp.cI(E)*sp.nu(I)*sp.tauSG*1e-3;
    mf.sGABA(I) = sp.cI(I)*sp.nu(I)*sp.tauSG*1e-3;
    mf.iGABA(E) = sp.g(rEG)*(mf.Vhld(E)-sp.VI)*mf.sGABA(E);
    mf.iGABA(I) = sp.g(rIG)*(mf.Vhld(I)-sp.VI)*mf.sGABA(I);
    
    mf.iSYN = mf.ixAMPA + mf.iAMPA + mf.iNMDA + mf.iGABA;
    
    mf.dphidmu(E) = rspnsFuncDeriv(sp.trp(E), tau(E), a(E), b(E), sigma(E), k2(E));
    mf.dphidmu(I) = rspnsFuncDeriv(sp.trp(I), tau(I), a(I), b(I), sigma(I), k2(I));
    mf.AA = mf.dphidmu./(S.*sp.gm).*(mu.*S-(S-1).*(mf.Vmn-sp.VL)).*sp.gm./sp.nu;
end


    function [ F ] = rootVmn( Vmn )
        J = 1 +  sp.gamma*exp(-sp.beta*Vmn);
        rho2(E) = sp.beta*sp.g(rEN)*sp.cE(E)*sp.tauSN*(Vmn(E)-sp.VE)*(J(E)-1)/(sp.gm(E)*J(E)^2)*1e-3;
        rho2(I) = sp.beta*sp.g(rIN)*sp.cE(I)*sp.tauSN*(Vmn(I)-sp.VE)*(J(I)-1)/(sp.gm(I)*J(I)^2)*1e-3;
        
        rho1(E) = sp.g(rEN)*sp.cE(E)*sp.tauSN/(sp.gm(E)*J(E))*1e-3;
        rho1(I) = sp.g(rIN)*sp.cE(I)*sp.tauSN/(sp.gm(I)*J(I))*1e-3;
        
        S(E) = 1 + TX(E)*sp.nuX(E) + TA(E)*sp.nu(E) + (rho1(E)+rho2(E))*sp.nu(E) + TG(E)*sp.nu(I);
        S(I) = 1 + TX(I)*sp.nuX(I) + TA(I)*sp.nu(E) + (rho1(I)+rho2(I))*sp.nu(E) + TG(I)*sp.nu(I);
        
        tau = sp.tau./S;
        mu(E) = (TX(E)*sp.nuX(E)+TA(E)*sp.nu(E)+rho1(E)*sp.nu(E))*(sp.VE-sp.VL)/S(E) + ...
            (rho2(E)*sp.nu(E)*(Vmn(E)-sp.VL)+TG(E)*sp.nu(I)*(sp.VI-sp.VL))/S(E);
        mu(I) = (TX(I)*sp.nuX(I)+TA(I)*sp.nu(E)+rho1(I)*sp.nu(E))*(sp.VE-sp.VL)/S(I) + ...
            (rho2(I)*sp.nu(E)*(Vmn(I)-sp.VL)+TG(I)*sp.nu(I)*(sp.VI-sp.VL))/S(I);
        
        F = mu + sp.VL - (sp.vThr-sp.vRst)*sp.nu.*tau*1e-3 - (mu+sp.VL-sp.vRst).*sp.trp.*sp.nu*1e-3 - Vmn;
    end
end