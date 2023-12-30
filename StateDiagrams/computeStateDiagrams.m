%   This code solves mean field (MF) and linear stability (LS) equations to
%   compute state diagrams. 
%
%   The code reproduces Fig.2A,B, the state diagrams shown in Fig.5A,B, 
%   and saves computed synaptic conductance parameters for the steady and
%   critical state primary networks in the "Parameters" subdirectory.

clear global
clear variables

global sp E I
global iNiG iThrE Vmn0 SYNC
global rEA rEG rEN xEA rIA rIG rIN xIA 

E = 1; I = 2;
[rEA, rEG, rEN, xEA, rIA, rIG, rIN, xIA] = deal(1,2,3,4,5,6,7,8);

% Set path to "Parameters" subdirectory
outFilePath = fullfile(pwd, 'StateDiagrams', 'Parameters');


%%%%%%%%%%%  Network Parameters  %%%%%%%%%%%

% Number of synaptic inputs per cell
sp.cX = [800, 800];
sp.cE = [800, 800];
sp.cI = [200, 200];

% Synaptic time constants
sp.tauD(E) = 1.0;   % ms, synaptic latency/delay
sp.tauD(I) = 1.0;

sp.tauAr = 0.2;     % ms, rising time constant
sp.tauNr = 2;
sp.tauGr = 0.5;

sp.tauAd = 2;       % ms, decay time constant
sp.tauNd = 100;
sp.tauGd = 5.0;

% Synaptic reversal potentials
sp.VE =  0;         % mV 
sp.VI = -70;

% Membrane time constants
sp.tau(E) = 20;      % ms
sp.tau(I) = 10;

sp.trp(E) = 2;       % ms
sp.trp(I) = 1;
                               
% Membrane potential constants
sp.vThr = -50;      % mV, firing threshold
sp.vRst = -55;      % mV, reset potential
sp.VL  = -70;       % mV, leakage potential

% Membrane capacitance and leakage conductance
sp.Cm(E) = 0.5;             % nF
sp.Cm(I) = 0.2;
sp.gm = sp.Cm./sp.tau*1e3;  % nS
    
% Voltage/magnesium NMDAR conductance dependence parameters
sp.beta = 0.062;       % 1/mV
sp.gamma = 1.00/3.57;  % mM/mM

% Normalization of time integral of gating variables
sp.tauSA = sp.tau(E);
sp.tauSN = sp.tau(E);
sp.tauSG = sp.tau(E);

% Population and external firing rates
nu0 = [5, 20];          % Hz
nuX0 = [nu0(E), nu0(E)];

% Current needed for E neuron to reach firing threshold
iThrE = -sp.gm(E)*(sp.vThr-sp.VL);

% Starting value for mean potential
V0 = (sp.vThr + sp.vRst)/2;
Vmn0 = [V0 V0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shft = 0;
ii = 1;

%
% Fig.2A: Compute state diagram in [iAMPA/iGABA - iX/iThr] plane
%
SYNC = 0;
sp.nu = nu0;
sp.nuX = nuX0;
x0ini = [0.0  1.2  0.5  0.4  0.0  1.0  0.4  0.3];
xx0ini = [0.3 0.2];
iNiG = 0.15;                % iNMDA/iGABA
sx = (0.56:-0.04:0.0).';    % iAMPA/iGABA
sy = (2.50:-0.08:0.98).';   % ixAMPA/iThrE

lambda = solve_gALL_lambda_fntwrk_iX(sx, sy, x0ini, xx0ini);

colorMap = jet;
figure(shft+ii);clf;set(gcf,'windowstyle','docked')
[~, hContour] = contourf(sx, sy, lambda.');
hContour.LineColor = 'none';
hContour.Fill = 'on';
hContour.LevelStep = 1;
colormap(colorMap)
caxis([-700, +300])
c = colorbar;
c.Label.String = 'lambda, 1/s';
ax=gca;
ax.TickDir = 'out';
ax.PlotBoxAspectRatio = [4.5 3 1];
xlabel('iAMPA / iGABA')
ylabel('iX / iThr')
grid on
hold on
[~, hContour] = contour(sx, sy, lambda.');
hContour.LineWidth = 1.5;
hContour.LineColor = 'w';   %
hContour.LevelList = 0.0;
ii = ii+1;

%
% Fig.2B: Compute oscillation frequency on the critical line 
%
SYNC = 1;
sp.nu = nu0;
sp.nuX = nuX0;
x0ini = [0.0  0.7  0.3  0.3  0.0  0.6  0.2  0.2  1.2];
xx0ini = 0;
sx = (0.0:0.01:0.56).';     % iAMPA/iGABA
sy = 0.15;                  % iNMDA/iGABA

[~, fntwrk] = solve_gALL_lambda_fntwrk_iX(sx, sy, x0ini, xx0ini);

figure(shft+ii);clf;set(gcf,'windowstyle','docked')
plot(sx, fntwrk, 'k', 'LineWidth', 1.5);
ax=gca;
ax.TickDir = 'out';
ax.PlotBoxAspectRatio = [4.5 3 1];
xlim([0, 0.56])
ylim([0, 200])
xlabel('iAMPA / iGABA')
ylabel('fntwrk, Hz')
ii = ii+1;

%
% Fig.5A,B: Compute state diagrams in [vX/vX* - gNMDA/gNMDA*] plane
%
xvar = (1.1:-0.005:0.9).';      % vX/vX*
yvar = (0:0.02:2).';            % gNMDA/gNMDA*
lambda = cell(2,1);

SYNC = 0;
sp.nu = nu0;
sp.nuX = nuX0;
x0ini = [0.0  1.2  0.5  0.4  0.0  1.0  0.4  0.3];
xx0ini = [0.3 0.2];
iNiG = 0.15;    % iN/iG=0.15, steady state primary network
sx = 0.2;       % iA/iG=0.20, steady state primary network
sy = 1.089;     % ixAMPA/iThrE=1.089, steady state primary network
[~, fntwrk, gSteadyAsync] = solve_gALL_lambda_fntwrk_iX(sx, sy, x0ini, xx0ini);
sp.g = round(sp.g,6);
spMF = sp;

% Save steady state primary network conductance parameters
fileName = 'paramSteadyNetwrk';
filePath = fullfile(outFilePath, fileName);
save(filePath, 'spMF');

SYNC = 0;
g0 = gSteadyAsync;
omega0 = 2*pi*fntwrk*1e-3;
lambda{1} = solve_gNMDA_lambda_fntwrk(xvar, yvar, SYNC, nuX0, g0, omega0, xx0ini);


SYNC = 1;
sp.nu = nu0;
sp.nuX = nuX0;
x0ini = [0.0  0.7  0.3  0.3  0.0  0.6  0.2  0.2  1.2];
xx0ini = 0;
sx = 0.4;       % iAMPA/iGABA
sy = 0.15;      % iNMDA/iGABA
[~, fntwrk, gOscilOnset, ixAMPAmf] = solve_gALL_lambda_fntwrk_iX(sx, sy, x0ini, xx0ini);
sp.g = round(sp.g,6);
spMF = sp;

% Save critical state primary network conductance parameters
fileName = 'paramCritclNetwrk';
filePath = fullfile(outFilePath, fileName);
save(filePath, 'spMF');

SYNC = 0;
g0 = gOscilOnset;
omega0 = 2*pi*fntwrk*1e-3;
xx0ini = [0; omega0];
lambda{2} = solve_gNMDA_lambda_fntwrk(xvar, yvar, SYNC, nuX0, g0, omega0, xx0ini);

colorMap = jet;
for i=1:2
    figure(shft+ii);clf;set(gcf,'windowstyle','docked')
    [~, hContour] = contourf(xvar, yvar, lambda{i}.');
    hContour.LineColor = 'none';
    hContour.Fill = 'on';
    hContour.LevelStep = 1;
    colormap(colorMap)
    c = colorbar;
    c.Label.String = 'lambda, 1/s';
    ax=gca;
    ax.TickDir = 'out';
    ax.PlotBoxAspectRatio = [4.5 3 1];
    xlabel('vX / vX*')
    ylabel('gNMDA / gNMDA*')
    grid on
    hold on
    [~, hContour] = contour(xvar, yvar, lambda{i}.');
    hContour.LineWidth = 2.5;
    hContour.LineColor = 'w';
    hContour.LevelList = 0.0;
    ii=ii+1;
end