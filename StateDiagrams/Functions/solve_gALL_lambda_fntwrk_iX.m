function [lambda, fntwrk, gmf, ixAMPAmf] = solve_gALL_lambda_fntwrk_iX(sx, sy, x0ini, xx0ini)
%   Solves MF & LS equations for synaptic conductances, instability growth
%   rate, and oscillation frequency given external and population firing
%   rates.
%
%   External and population rates are passed through 
%   global variables sp.nuX and sp.nu, repsectively
%
%   Case SYNC = 0, general solution
%       Input
%           sx      - AMPA / GABA current balance
%           sy      - external input / firing threshold current balance
%           x0ini   - starting values for 8 synaptic conductances 
%           xx0ini  - starting values for instability growth rate and oscillation frequency 
%       Output
%           gmf     - 8 synaptic conductances
%           fntwrk  - frequency of oscillation
%           lambda  - oscillatory instability growth rate 
%           ixAMPAmf- external AMPA current
%
%   Case SYNC = 1, solution for critical point lambda = 0
%       Input
%           sx      - AMPA / GABA current balance
%           sy      - NMDA / GABA current balance
%           x0ini   - starting values for 8 synaptic conductances and oscillation frequency 
%       Output
%           gmf     - 8 synaptic conductances
%           fntwrk  - frequency of oscillation
%           ixAMPAmf- external AMPA current

global sp mf E
global iAiG iNiG
global iThrE iXiT SYNC

if SYNC
    neq = 9;
else
    neq = 8;
end
Nx = length(sx);
Ny = length(sy);
gmf = zeros(Nx, Ny,8);
fntwrk = zeros(Nx, Ny);
lambda = zeros(Nx, Ny);
ixAMPAmf = zeros(Nx, Ny,2);
x0 = zeros(Nx,neq);
xx0 = zeros(Nx,2);
x0(1,:) = x0ini;
xx0(1,:) = xx0ini;

options = optimoptions(@fsolve,'Display','off');
for j = 1:Ny
    for i = 1:Nx
        iAiG = sx(i);
        if SYNC
            iNiG = sy(j);
        else
            iXiT = sy(j);
        end
        [x, ~, exitflag] = fsolve(@root_gALL_omega, x0(i,:), options);
        if exitflag < 1
            fprintf('No solution for g and omega, i=%d, Exitflag = %d\n', i, exitflag)
            return
        end
        x = abs(x);
        % Unpack x
        if SYNC
            g = x(1:end-1);
            omega = x(end)*1e3;
            fntwrk(i,j) = omega/(2*pi);    % Hz
        else
            g = x;
        end
        % Save current solution as initial condition for the next iteration cycle
        x0(i,:) = x;
        if i <= Nx
            x0(i+1,:) = x;
        end
        [~, mf] = mfa(sp);
        if ~SYNC
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
                if i <= Nx
                    xx0(i+1,:) = xx;
                end
            end
        end
        gmf(i,j,:) = g;
        ixAMPAmf(i,j,:) = mf.ixAMPA;
    end
    if ~SYNC
        fprintf('iX/iT=%.2f,  lambda=%.1f,  fntwrk=%.1f Hz\n', sy(j), lambda(i,j), fntwrk(i,j)');
    else
        fprintf('iN/iG=%.2f,  iXE/iTE=%.3f,  fntwrk=%.1f Hz\n', sy(j), mf.ixAMPA(E)/iThrE, fntwrk(i,j)');
    end
end
end