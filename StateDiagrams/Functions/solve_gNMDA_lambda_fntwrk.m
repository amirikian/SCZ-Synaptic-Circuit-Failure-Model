function [lambda, fntwrk, gN] = solve_gNMDA_lambda_fntwrk(xvar, yvar, sync, nuX0, g0, omega0, xx0ini)
%   Solves MF & LS equations for population firing rates, instability growth rate, 
%   and oscillation frequency given conductances and external firing rate.
%
%   Population rates are passed through global variable sp.nu
%   
%   Case sync = 0, general solution
%       Input
%           xvar  - modulation strength of external rates
%           yvar  - modulation strength of NMDAR conductances
%           xx0ini- starting values for lambda and oscillation frequency 
%       Output
%           lambda      - oscillatory instability growth rate 
%           fntwrk      - frequency of oscillation
%           gN          - NMDAR conductances
% 
%   Case sync = 1, solution for critical point lambda = 0
%       Input
%           xvar  - modulation strength of external rates
%           nuX0  - base external rates
%           g0    - starting values for NMDAR condactances
%           omega0- starting value for oscillation frequency 
%       Output
%           gN    - NMDAR conductances
%           fntwrk- frequency of oscillation

global sp mf E I
global rEN rIN
global grEN0 grIN0

nx = length(xvar);
ny = length(yvar);
vx = zeros(nx,ny,2);
gN = zeros(nx,ny,2);
fntwrk = zeros(nx,ny);
lambda = zeros(nx,ny);
gNnuOmega0 = [1, sp.nu, omega0];
nu0 = repmat(sp.nu,ny,1);
xx0 = zeros(nx,2);
xx0(1,:) = xx0ini;

sp.g = g0;
grEN0 = g0(rEN);
grIN0 = g0(rIN);

options = optimoptions(@fsolve,'Display','off');
for i = 1:nx
    for j = 1:ny
        vx(i,j,:) = xvar(i)*nuX0;
        sp.nuX = vx(i,j,:);
        if sync
            [gNnuOmega, ~, exitflag] = fsolve(@root_gNMDA_nu_omega, gNnuOmega0, options);
            gN(i,j,:) = gNnuOmega(1)*[grEN0,grIN0];
            v = gNnuOmega(2:3);
            omega = gNnuOmega(end)*1e3;
            fntwrk(i) = omega/(2*pi);
            gNnuOmega0 = gNnuOmega;
        else
            sp.g(rEN) = yvar(j)*grEN0;
            sp.g(rIN) = yvar(j)*grIN0;
            gN(i,j,:) = yvar(j)*[grEN0,grIN0];
            [v, ~, exitflag] = fsolve(@root_nu, nu0(j,:), options);
        end
        if exitflag < 1
            [v, ~, exitflag] = fsolve(@root_nu, nu0(j+1,:), options);
            if exitflag < 1
                fprintf('No solution, sync = %d, Exitflag = %d\n', sync, exitflag)
                return
            else
                nu0(j,:) = v;
                if j <= ny
                    nu0(j+1,:) = v;
                end
            end
        else
            nu0(j,:) = v;
            if j <= ny
                nu0(j+1,:) = v;
            end
        end
        [~, mf] = mfa(sp);
        if ~sync
            [xx, ~, exitflag] = fsolve(@root_lambda_omega, xx0(i,:), options);
            if exitflag < 1
                fprintf('No solution for lambda and omega, Exitflag = %d\n', exitflag)
                return
            else
                % Unpack xx
                omega = abs(xx(2)*1e3);
                fntwrk(i,j) = omega/(2*pi);
                lambda(i,j) = xx(1)*1e3;
                % Save current solution as initial condition for the next iteration cycle
                xx0(i,:) = xx;
                if i <= nx
                    xx0(i+1,:) = xx;
                end
            end
        end
    end
    fprintf('vX/vX*=%.3f,  vE=%.1f,  vI=%.1f,  fntwrk=%.1f Hz,  lambda=%.1f\n', xvar(i), v(E), v(I), fntwrk(i,j), lambda(i,j));
end
end