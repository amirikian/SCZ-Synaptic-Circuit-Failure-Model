function dphidmu = rspnsFuncDeriv(trp, tau, a, b, sigma, k2)
%   Derivative of neuron's current-frequency response function
%   
%     Input:
%         trp    - refractory period time
%         tau    - membrane time constant
%         sigma  - standard deviation of synaptic noise
%         a,b,k2 - constants provided in the parent function mfa.m
%     Output:
%         dphidmu- derivative of response function

nu = rspnsFunc(trp, tau, a, b);
if nu
    dphidmu = nu*nu*tau*1e-3/sigma*sqrt(pi)*((1+0.5*k2)*erfcx(-a)-erfcx(-b));
else
    dphidmu = 0;
end
end