function i = intgrtDynEqs(odefun, i0)
%   Integrates dynamical equations for the membrane potentials  
%   and NMDAR gating variables until the membrane potential of a neuron 
%   crosses the firing threshold.
%
%       Input
%           odefun - r.h.s. of the dynamical equations for the membrane 
%                    potentials  and NMDAR gating variables
%           i0     - time step at which integration starts
%       Output
%           i      - time step at which integration stops 

global neq n hi 
global lastSpk trp
global tspan N Y
global sp

F = zeros(neq,2);

for i = i0+1:N
  ti = tspan(i-1);
  yi = Y(:,i-1);
  F(:,1) = odefun(ti,yi);
  F(:,2) = odefun(ti+hi,yi+hi*F(:,1));
  Y(:,i) = yi + (hi/2)*(F(:,1) + F(:,2));

% Set Y of neurons in refractory period to vRst
  k = ti+hi <= trp+lastSpk;
  if any(k)
      Y(k,i) = sp.vRst;
  end
  
  k = (ti+hi > trp+lastSpk) & (ti <= trp+lastSpk);
  if any(k)
      trpend = lastSpk(k) + trp(k);
      Y(k,i) = sp.vRst + (ti+hi-trpend).*(F(k,1)+F(k,2))/2;
  end

  if(any(Y(1:n,i) > sp.vThr))
      break;
  end
end