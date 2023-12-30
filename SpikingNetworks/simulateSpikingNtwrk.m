%   This code simulates the dynamics of recurrent spiking network model for 
%   specified AMPAR, NMDAR, and GABAR conductances and saves the spike 
%   trains of individual neurons in the "Simulations" subdirectory.
%
%   The values of the synaptic conductances can be either loaded from the 
%   "Parameters" subdirectory using one of the files created previously by 
%   the computeStateDiagrams.m, or specified directly by using/modifying 
%   conductance parameters for the steady and critical state primary networks 
%   listed in this code. 
%
%   To simulate networks corresponding to Fig.3A1, Fig.3A2, Fig.3B1 or 
%   Fig.3B2, first run the computeStateDiagrams.m to compute and save 
%   conductances corresponding to the steady and critical state primary 
%   networks. Then, run this code with the ntwrk2simul variable set to 
%   'Fig.3A1', 'Fig.3A2', 'Fig.3B1', or 'Fig.3B2'.
%
%   Alternatively, set ntwrk2simul='myNtwrk' and use the values of synaptic
%   conductances for the steady or critical state primary networks
%   explicitly specified in this code. Behavior of the network can be
%   explored by varying the modulation strengths of external inputs fvX 
%   and/or NMDAR conductances fgN.

clear global
clear variables

% Initialize random number generator
iseed = 12345;
rng(iseed);
rngCur = rng;

% Set path to "Parameters" subdirectory
inpFilePath = fullfile(pwd, 'StateDiagrams', 'Parameters');

% Allocate global variables
global wEE wEI wIE wII
global vEindx vIindx sNindx
global n neq hi tspan trp N Y
global lastSpk dlySpkTimesE dlySpkTimesI Sext
global rEA rEG rEN xEA rIA rIG rIN xIA E I sp

E = 1; I = 2;
[rEA, rEG, rEN, xEA, rIA, rIG, rIN, xIA] = deal(1,2,3,4,5,6,7,8);

% Choose network to simulate:
% ntwrk2simul = 'Fig3A1' - network corresponding to Fig.3A1 in the paper
% ntwrk2simul = 'Fig3B1' - network corresponding to Fig.3B1 in the paper
% ntwrk2simul = 'Fig3A2' - network corresponding to Fig.3A2 in the paper
% ntwrk2simul = 'Fig3B2' - network corresponding to Fig.3B2 in the paper
% ntwrk2simul = 'myNtwrk' - network with parameters specified below

% ntwrk2simul = 'Fig3A1';
% ntwrk2simul = 'Fig3B1';
% ntwrk2simul = 'Fig3A2';
ntwrk2simul = 'Fig3B2';
% ntwrk2simul = 'myNtwrk';

switch ntwrk2simul
    case 'Fig3A1'
        load(fullfile(inpFilePath, 'paramSteadyNetwrk'), 'spMF');
        load(fullfile(inpFilePath, 'rnginiSteadyNetwrk_Fig3A1'), 'rngCur');
        outFileName = 'simulSteadyNetwrk_';
        sp = spMF;
        sp.fgN = 1.0;   % Modulation strength of NMDAR conductances
        sp.fvX = 1.0;   % Modulation strength of external spike rate
    case 'Fig3B1'
        load(fullfile(inpFilePath, 'paramSteadyNetwrk'), 'spMF');
        load(fullfile(inpFilePath, 'rnginiSteadyNetwrk_Fig3B1'), 'rngCur');
        outFileName = 'simulSteadyNetwrk_';
        sp = spMF;
        sp.fgN = 1.0;    % Modulation strength of NMDAR conductances
        sp.fvX = 1.05;   % Modulation strength of external spike rate
    case 'Fig3A2'
        load(fullfile(inpFilePath, 'paramCritclNetwrk'), 'spMF');
        load(fullfile(inpFilePath, 'rnginiCritclNetwrk_Fig3A2'), 'rngCur');
        outFileName = 'simulCritclNetwrk_';
        sp = spMF;
        sp.fgN = 1.0;   % Modulation strength of NMDAR conductances
        sp.fvX = 1.0;   % Modulation strength of external spike rate
    case 'Fig3B2'
        load(fullfile(inpFilePath, 'paramCritclNetwrk'), 'spMF');
        load(fullfile(inpFilePath, 'rnginiCritclNetwrk_Fig3B2'), 'rngCur');
        outFileName = 'simulCritclNetwrk_';
        sp = spMF;
        sp.fgN = 1.0;    % Modulation strength of NMDAR conductances
        sp.fvX = 1.05;   % Modulation strength of external spike rate
    case 'myNtwrk'
        %%%%%%%%%%%  Specify Basic Network Parameters  %%%%%%%%%%%
        
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
        
        % Membrane capacitance
        sp.Cm(E) = 0.5;             % nF
        sp.Cm(I) = 0.2;
        
        % Voltage/magnesium NMDAR conductance dependence parameters
        sp.beta = 0.062;       % 1/mV
        sp.gamma = 1.00/3.57;  % mM/mM
       
        % Normalization of time integral of gating variables
        sp.tauSA = sp.tau(E);
        sp.tauSN = sp.tau(E);
        sp.tauSG = sp.tau(E);

        % Population spike rates
        sp.nu(E) = 5;   % Hz
        sp.nu(I) = 20;  % Hz
        
        % External input spike rates
        sp.nuX = [sp.nu(E) sp.nu(E)];
                       
        % Default modulation parameters
        sp.fgN = 1;     % Modulation strength of NMDAR conductances
        sp.fvX = 1;     % Modulation strength of external spike rate
        
        %%%%%%%%%%%  Specify Synaptic Conductnaces  %%%%%%%%%%%
        
        %%% Use critical state primary network conductance parameters
        %%%     iN/iG = 0.15, iA/iG = 0.4, iXE/iTE=1.089
        sp.g(rEA) = 0.019317;  % nS,     
        sp.g(rEG) = 0.143892;
        sp.g(rEN) = 0.059546;
        sp.g(xEA) = 0.129921;
        sp.g(rIA) = 0.015856;
        sp.g(rIG) = 0.118723;
        sp.g(rIN) = 0.049058;
        sp.g(xIA) = 0.106641;
        outFileName = 'simulCritcl_';
        
        % Use Fig3B2 rng settigs
        load(fullfile(inpFilePath, 'rnginiCritclNetwrk_Fig3B2'), 'rngCur');
        % Use Fig3B2 modulations
        sp.fgN = 1.00;  % Modulation strength of NMDAR conductances
        sp.fvX = 1.05;  % Modulation strength of external spike rate
        
%         % Use Fig3A2 rng settigs
%         load(fullfile(inpFilePath, 'rnginiCritclNetwrk_Fig3A2'), 'rngCur');
%         % Use Fig3A2 modulations
%         sp.fgN = 1.0;   % Modulation strength of NMDAR conductances
%         sp.fvX = 1.0;   % Modulation strength of external spike rate
        
%         %%% Use steady state primary network conductance parameters
%         %%%     iN/iG = 0.15, iA/iG = 0.2, iXE/iTE=1.089
%         sp.g(rEA) = 0.006722;      % nS
%         sp.g(rEG) = 0.100341;
%         sp.g(rEN) = 0.041501;
%         sp.g(xEA) = 0.129802;
%         sp.g(rIA) = 0.005513;
%         sp.g(rIG) = 0.082773;
%         sp.g(rIN) = 0.034178;
%         sp.g(xIA) = 0.106452;
%         outFileName = 'simulSteady_';
        
%         % Use Fig3B1 rng settigs
%         load(fullfile(inpFilePath, 'rnginiSteadyNetwrk_Fig3B1'), 'rngCur');
%         % Use Fig3B1 modulations
%         sp.fgN = 1.00;  % Modulation strength of NMDAR conductances
%         sp.fvX = 1.05;  % Modulation strength of external spike rate
        
%         % Use Fig3A1 rng settigs
%         load(fullfile(inpFilePath, 'rnginiSteadyNetwrk_Fig3A1'), 'rngCur');
%         % Use Fig3A1 modulations
%         sp.fgN = 1.00;  % Modulation strength of NMDAR conductances
%         sp.fvX = 1.00;  % Modulation strength of external spike rate
        
%         % Use my modulation parameters
%         sp.fgN = 1.00;  % Modulation strength of NMDAR conductances
%         sp.fvX = 1.05;  % Modulation strength of external spike rate
        
        % Label file name
        fileLabel = sprintf('fgN=%.2f_fvX=%.2f_', sp.fgN, sp.fvX);
        outFileName = strcat(outFileName, fileLabel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        error('Undefined network')
end
rng(rngCur);
sp.rng = rngCur;

% Output file name
outFileName = strcat(outFileName, ntwrk2simul, '.mat');

% Network size and composition
n = 5000;
sp.nE = round(0.8*n);
sp.nI = n-sp.nE;

% Connection probability
pcon = (sp.cE(E)+sp.cI(E))/n;

% Connectivity matrices
wEE = binornd(1, pcon, sp.nE, sp.nE);
wEI = binornd(1, pcon, sp.nE, sp.nI);
wIE = binornd(1, pcon, sp.nI, sp.nE);
wII = binornd(1, pcon, sp.nI, sp.nI);

% Adjust NMDAR conductances and external spike rates
sp.g(rEN) = sp.fgN*sp.g(rEN);
sp.g(rIN) = sp.fgN*sp.g(rIN);
sp.nuX = sp.fvX*sp.nuX;
rmnext = sp.cX(E)*[sp.nuX(E)*ones(sp.nE,1); sp.nuX(I)*ones(sp.nI,1)]*1e-3; %   spks/ms

% Allocate memory for membrane potentials and gating variables.
y0 = zeros(n+sp.nE, 1);
neq = length(y0);

% Initialize membrane potentials and gating variables
%  sy0 - randomization scale of membrane potentials, 0 <= sy0 < 1
sy0 = 0;    % in this simulation initial potentials are not randomized
y0(1:n) = 0.5*(sp.vThr+sp.vRst+sy0*(sp.vThr-sp.vRst)*(1-2*rand(1,n)));
y0(sNindx) = 0;

% Create storage to keep spike times
maxNumSpks = 1000;
allSpkTimes = -inf(maxNumSpks, n);

maxNumSpksE = 15;
maxNumSpksI = 35;
dlySpkTimesE = -inf(maxNumSpksE, sp.nE);
dlySpkTimesI = -inf(maxNumSpksI, sp.nI);

% Set simulation time and integration step
T = 3500;       %   ms
dt = 0.1;       %   ms
sp.T = T;
sp.tstp = dt;
hi = dt;
tspan = 0:dt:T;
N = length(tspan);

vEindx = 1:sp.nE;
vIindx = sp.nE+1:n;
sNindx = n+1:n+sp.nE;

sp.ensmE = (vEindx).';
sp.ensmI = (vIindx).';

% Generate Poisson spike trains and gating variables for external inputs

% Define simulation protocol
extEvnt = 'trialON';
tExtEvnt = 0;
tevnt = -100*sp.tauAd;
xSpkTimes = cell(n,1);
sp.tTrlON = tExtEvnt;
sp.tTrlOFF = T;

while tevnt < T
    if (tevnt >= tExtEvnt )
        switch (extEvnt)
            case 'trialON'
                tExtEvnt = sp.tTrlOFF;
        end
    end
    nXspks = poissrnd(rmnext*(tExtEvnt-tevnt), [n, 1]);
    for icell=1:n
        arrvlTimes = tevnt + (tExtEvnt-tevnt)*rand(nXspks(icell),1);
        if ~isempty(arrvlTimes)
            xSpkTimes{icell} = [xSpkTimes{icell}; sort(arrvlTimes)];
        end
    end
    tevnt = tExtEvnt;
end

% Parallel computing of gating variables for external currents
nN = n;
nT = round(T/dt);
sext = zeros(n, nT+1);

parfor icell = 1:nN
    for it = 1:nT+1
        t = (it-1)*dt;
        [idx1, idx2] = findInSorted(xSpkTimes{icell}, [t-10*sp.tauAd, t]);
        sext(icell, it) = sum(exp((xSpkTimes{icell}(idx1:idx2)-t)/sp.tauAd) ...
            - exp((xSpkTimes{icell}(idx1:idx2)-t)/sp.tauAr));
    end
end
Sext = sext*sp.tauSA/(sp.tauAd-sp.tauAr);

% Generate spike times of the preceding t=0 external input spikes
nXspks = [ poissrnd(10*max(sp.tauNr,sp.tauAd)*sp.nu(E)*1e-3, [sp.nE, 1]);...
           poissrnd(10*sp.tauGd*sp.nu(I)*1e-3, [sp.nI, 1]) ];
for icell = vEindx
    arrvlTimes = -10*max(sp.tauNr,sp.tauAd)*sort(rand(nXspks(icell),1));
    dlySpkTimesE(1:length(arrvlTimes), icell) = arrvlTimes;
end
for icell = vIindx
    arrvlTimes = -10*sp.tauGd*sort(rand(nXspks(icell),1));
    dlySpkTimesI(1:length(arrvlTimes), icell-sp.nE) = arrvlTimes;
end


% Set integration options
it = 1;
t0 = tspan(it);
lastSpk = [dlySpkTimesE(1,:), dlySpkTimesI(1,:)].';
tDsp = 0;
tStp = 250; % ms, time interval for displying the number of spiking cells

spkForgetTime = zeros(1, n);
spkForgetTime(vEindx) = max(10*sp.tauAd, 10*sp.tauNr);
spkForgetTime(vIindx) = 10*sp.tauGd;

tau = zeros(n,1);
trp = zeros(n,1);
tau(vEindx) = sp.tau(E);
tau(vIindx) = sp.tau(I);
trp(vEindx) = sp.trp(E);
trp(vIindx) = sp.trp(I);

tD = repelem(sp.tauD, [sp.nE, sp.nI]);

Y = zeros(neq,N);
Y(:,1) = y0;
tic;
while t0 < T
    nt = intgrtDynEqs(@rhsDynEqs, it);
    
    it = nt;
    
    % Find neurons that crossed firing threshold and save their spiking times
    spkCells = find(Y(1:n,nt)>sp.vThr);
    te = (tspan(it-1) + dt*(sp.vThr - Y(spkCells,nt-1))./(Y(spkCells,nt) - Y(spkCells,nt-1))).';
    
    lastSpk(spkCells) = te;
    
    if all(isinf(allSpkTimes(end, spkCells)))
        allSpkTimes(end, spkCells) = te;
        allSpkTimes(:, spkCells) = circshift(allSpkTimes(:, spkCells), [1 0]);
    else
        fprintf('Some earlier spikes are lost');
        allSpkTimes(end, spkCells) = te;
        allSpkTimes(:, spkCells) = circshift(allSpkTimes(:, spkCells), [1 0]);
    end
    
    spkCellsE = spkCells(spkCells<=sp.nE);
    if ~isempty(spkCellsE)
        teE = te(1:length(spkCellsE));
    a = teE + tD(spkCellsE) - dlySpkTimesE(end, spkCellsE);
    if any(a < spkForgetTime(spkCellsE))
        error('Too small buffer for storing delayed E spike times');
    end
    dlySpkTimesE(end, spkCellsE) = teE + tD(spkCellsE);
    dlySpkTimesE(:, spkCellsE) = circshift(dlySpkTimesE(:, spkCellsE), [1 0]);
    end
    spkCellsI = spkCells(spkCells>sp.nE)-sp.nE;
    if ~isempty(spkCellsI)
        teI = te(end-length(spkCellsI)+1:end);
    a = teI + tD(spkCellsI+sp.nE) - dlySpkTimesI(end, spkCellsI);
    if any(a < spkForgetTime(spkCellsI+sp.nE))
        error('Too small buffer for storing delayed I spike times');
    end
    dlySpkTimesI(end, spkCellsI) = teI + tD(spkCellsI+sp.nE);
    dlySpkTimesI(:, spkCellsI) = circshift(dlySpkTimesI(:, spkCellsI), [1 0]);
    end
    
    % Display the number of spiking cells in the present threshold crossing event
    if (min(te) >= tDsp)
        fprintf('t = %.4f,  nSpikingCells = %d\n', min(te), length(spkCells));
        tDsp = tDsp + tStp;
    end
         
    % Reset membrane potential
    Y(spkCells,nt) = sp.vRst;
    t0 = tspan(it);
end
toc;

% Save simulation results
outFilePath = fullfile(pwd, 'SpikingNetworks', 'Simulations');
filePath = fullfile(outFilePath, outFileName);
save(filePath, 'sp', 'allSpkTimes', '-v7.3');