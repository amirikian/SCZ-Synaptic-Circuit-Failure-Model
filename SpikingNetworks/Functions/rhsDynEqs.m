function ydot = rhsDynEqs( t, y )
%   Computes r.h.s. of the dynamical equations (rate of change) for the 
%   membrane potentials and NMDAR gating vriables.
%
%       Input
%           t    - time
%           y    - membrane potentials and NMDAR gating variables at time t
%       Output
%           ydot - rate of change of the membrane potentials and NMDAR 
%                  gating variables at time t

global wEE wEI wIE wII
global dlySpkTimesE dlySpkTimesI Sext
global vEindx vIindx sNindx
global rEA rEG rEN xEA rIA rIG rIN xIA 
global sp E I hi

% Unpack variables
vE = y(vEindx);
vI = y(vIindx);
sN = y(sNindx);

% Calculate currents and potentials
sX = Sext(:,1+round(t/hi));
iEX = sp.g(xEA)*(vE-sp.VE).*sX(vEindx);
iIX = sp.g(xIA)*(vI-sp.VE).*sX(vIindx);

arrvdSpkTimes = t - dlySpkTimesI;
arrvdSpkTimes(arrvdSpkTimes<0) = inf;
sG = sp.tauSG/(sp.tauGd-sp.tauGr)*(sum(exp(-arrvdSpkTimes/sp.tauGd)-exp(-arrvdSpkTimes/sp.tauGr))).';
arrvdSpkTimes = t - dlySpkTimesE;
arrvdSpkTimes(arrvdSpkTimes<0) = inf;
sA = sp.tauSA/(sp.tauAd-sp.tauAr)*(sum(exp(-arrvdSpkTimes/sp.tauAd)-exp(-arrvdSpkTimes/sp.tauAr))).';
xN = sp.tauSN/sp.tauNr*sum(exp(-arrvdSpkTimes(:,vEindx)/sp.tauNr)).';

iEG = sp.g(rEG)*(vE-sp.VI).*(wEI*sG);       %   pA
iIG = sp.g(rIG)*(vI-sp.VI).*(wII*sG);

iEN = sp.g(rEN)*(vE-sp.VE)./(1+sp.gamma*exp(-sp.beta*vE)).*(wEE*sN);
iIN = sp.g(rIN)*(vI-sp.VE)./(1+sp.gamma*exp(-sp.beta*vI)).*(wIE*sN);

iEA = sp.g(rEA)*(vE-sp.VE).*(wEE*sA);
iIA = sp.g(rIA)*(vI-sp.VE).*(wIE*sA);

iampa = [iEA; iIA];
inmda = [iEN; iIN];
igaba = [iEG; iIG];
ixampa= [iEX; iIX];

itot = iampa + inmda + igaba + ixampa;

vEdot = -(vE-sp.VL)/sp.tau(E) - itot(vEindx)/sp.Cm(E)*1e-3;
vIdot = -(vI-sp.VL)/sp.tau(I) - itot(vIindx)/sp.Cm(I)*1e-3;
sNdot = (-sN + xN)/sp.tauNd;

% Pack variables
ydot = [vEdot; vIdot; sNdot];

end