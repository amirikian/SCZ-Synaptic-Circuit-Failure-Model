function nu = rspnsFunc(trp, tau, a, b)
%   Neuron's current-frequency response function
%   
%     Input:
%         trp - refractory period time
%         tau - membrane time constant
%         a,b - limits of integration provided in the parent function mfa.m
%     Output:
%         nu  - firing rate

nu = 1e3/(trp + sqrt(pi)*tau*integral(@(x) erfcx(-x), b, a, 'ArrayValued', 1));
end