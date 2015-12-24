function [scenario] = simulationParams_eh(inputParams)
% simulationParams_eh.m
%
% USAGE:
%       simulationParam_eh(inputParams)
%
%  inputParams - vector with the following optional subfields:
%                optimizerFnc
%                objType
%                oracleMode
%                numSensors
%                Horizon
%                nodeEnergyInitScale
%                BW
%                dopplerFreqMax
%
% defaults exist, so these are not necessary
% 
% parameter setup for WSN with energy harvesting simulations
%
%
% Copyright (C) 2011  Nick Roseveare
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% By: Nick Roseveare, 10-2011
global sensorData recycleMCdraws


%% ================ set defaults if no inputs ================
optimizerFnc = @optimizeQuant_std;
objType    = 'BLUEwfwf';% can be sumrate or BLUEwfwf or BLUErwfrwf (water-filling vs. reverse-WF)
oracleMode = 0;
singleShot = 0;
stopIferror= 1;
numSensors = 5;
nsamples   = 10;
Horizon    = 3;
BW         = 60;
newWSNposEachMCrun = 0;
newTrajEachMCrun   = 0;
dopplerFreqMax = 0.1; % good values for tstep=1.5 are: cyclic corr=>0.1, uncorr=>10, highly corr=> 0.01
%% ======= read in input values (which can override the defaults) =========
if nargin
    strs = fields(inputParams);
    for ss = 1:length(strs);
        if any(strcmp(strs{ss},{'optimizerFnc','objType','oracleMode','stopIferror','numSensors','nsamples','Horizon','nodeEnergyInitScale','BW','dopplerFreqMax','singleShot','newWSNposEachMCrun','newTrajEachMCrun'}))
            eval(sprintf('%s = inputParams.%s;',strs{ss},strs{ss}));
        else
            error('Input subfield \''%s\'' is invalid\n',strs{ss});
        end
    end
end
% set the default scaling of initial battery 
if nargin == 0
    nodeEnergyInitScale = ones(numSensors,1);
else
    strs = fields(inputParams);
    if any(strcmp(strs{ss},'nodeEnergyInitScale'))
        nodeEnergyInitScale = ones(numSensors,1);
    else
        nodeEnergyInitScale = ones(numSensors,1); 
    end
end

%% ============== Control whether debug stops occur on errors =============
if stopIferror
    dbstop if error 
else 
    dbclear if error
end

%% ================= Simulation Setup parameters =========================
scenario.descrp = 'This structure contains most of the settings for the generation of the data, sensor and communication parameters used for the simulation';
scenario.info   = 'There are other variables saved in the same data file, but they do not lend themselves\nto easy adding of text lines: "reportsAtFC" contains all of the results and data\nused to generate the plots, and sensorState contains a few bits of data not here.';
scenario.date   = datestr(date,'yyyy-mm-dd'); % just for reference in case the date is not recorded elsewhere
% Data and plotting controls
scenario.titleBars        = 1;
scenario.saveDataAndPlots = 0; % indicate whether all data will be saved
scenario.makePlots        = 1; % data will only be saved if above true, this determines if plots are shown; if not shown, they will not be saved
scenario.checkForOW       = 0; % check for overwriting previous files
scenario.plotExtraSetupInfo = 1;
scenario.plotOptmzProgress = 0;
scenario.saveOptmzProgress = 0;
scenario.plotOptmzProgOnline = 0;

% Simulation and debug controls
scenario.oracleMode       = oracleMode;
scenario.stopiferror      = stopIferror;
scenario.singleShot       = singleShot;
scenario.useLastTruth     = 1;
scenario.useLastWSNpos    = 1;
scenario.useFirstNew      = 1;
scenario.newTrajEachMCrun = newTrajEachMCrun;
scenario.newWSNposEachMCrun = newWSNposEachMCrun;
scenario.simDataSaveDirectory = [pwd, '/'];%'../../PhdWork/Code/';%'./'; %
scenario.randDistWSNpos   = 1;

% Objectives and metrics
scenario.objectiveType = objType;
scenario.lifetimeMetric = 'power';% can be 'power' or 'mse'
scenario.vectorMetricType = 'onedimen';% options: 'trace','vector','onedimen'
scenario.useDimen = 1; % 1 = x position, 2 = y position, 3 = x velocity, etc.

% Scenario parameters
scenario.fixNoiseSamples = 0;
scenario.tsample  = 1;
scenario.nsamples = nsamples;
scenario.numSamplesPredBased = 3; % for non-single-step Markov processes
scenario.numDimen = 2;
scenario.observRegion = [0 20 0 20 7];% xmin xmax ymin ymax vmagmax
scenario.runIdx = 1; % default mc run index for non-MC runs to store noise values when fixNoiseSamples = True
% scenario.lifetimeMetric   = 'power';% can be 'power' or 'mse'
% scenario.alpha            = 1; % WARNING: scenario may change this; if change made here, make sure that later scenario.alpha set locations are also changed
% scenario.plotIndivNode2D = 0;
% scenario.mseLevel = 8;
% scenario.reuseScenarioSetup = 0; % WARNING: if 1 => above should be 1: i.e. if this is '1' it can be overridden by fixNoiseSamples if it is '0'
% scenario.lifeTimeThres    = 0.75;


recycleMCdraws.defaultRunIdx = scenario.runIdx;
if isfield(recycleMCdraws,'reuseNoise')
    scenario.fixNoiseSamples = recycleMCdraws.reuseNoise;
else
    recycleMCdraws.reuseNoise = scenario.fixNoiseSamples;
end

time = 0:scenario.tsample:(scenario.tsample*(scenario.nsamples-1));
scenario.time = time;

dim     = scenario.numDimen;
numElem = dim*2; % number of elements in the state vector PV * dim
scenario.numElem = numElem;

if scenario.singleShot
     Horizon = scenario.nsamples;
end

% Find dynamic range of truth: multiply by .95 to give better resolution to
% points on the interior (there are fewer points near the outside
% (typically); multiply by 1/2 because W represents the range [-W,W]
Wpos = 0.95*(0.5)*max([scenario.observRegion(2)-scenario.observRegion(1)  scenario.observRegion(4)-scenario.observRegion(3)]);
Wvel = scenario.observRegion(5);

%% Individual sensor node parameters
scenario.numSensors = numSensors;
scenario.useRealisticRngAng = 1;
scenario.noiseVar = 5*ones(numSensors,1); %(normally 1, smaller for "nice" plots) base measurement noise variance
scenario.angNoiseVar = .05*pi;
scenario.rangeNoiseVar = .1; % x distance
scenario.q          = 1.5;
scenario.p0         = 40;%std->had on 4 % scale the initial covariance matrix by this no. (needs to start large to be reasonable)
scenario.scaleFunction = @reciprocalUnitModulus; % the heuristic "resource policy" scaling function for the power-bit usage control
scenario.FCnodePosition = [mean(scenario.observRegion(1:2));mean(scenario.observRegion(3:4))];

% Generate the true state dynamics
scenario.qtrue      = 1;
scenario.initState  = [0*ones(1,scenario.numDimen) (2+.5*Wpos/scenario.nsamples)*ones(1,scenario.numDimen)]';
disp('Generating system trajectory')
scenario.trueSys = setTrueTrajectory(scenario);

% Get sensor positions and setup the sensor KF states
disp('Generating sensor network setup')
sensorState = makeSensorStructEH(scenario,time);
scenario.sensorState = sensorState;
scenario.sensorPositions = [sensorState.senPos];

[F,Q] = makeFQ(dim,scenario.tsample,scenario.qtrue);
H = [eye(dim) zeros(dim)];
Pnext = F*sensorState(1).P*F' + Q;
scenario.nominalP = diag(Pnext - (Pnext*H'/(scenario.noiseVar(1)*eye(dim)+H*Pnext*H'))*H*Pnext);


%% Comm system characteristics 
commParam.n0         = 0.5*ones(1,numSensors); %std 1*...: %1*ones(1,numSensors);% noise floor
commParam.PrBitError = 0.05; % probability of bit error
commParam.BW         = BW;%std=60; % 65 % total bandwidth allowed to be used in a single reporting period
commParam.BWmax      = 25; % so that optimization results don't mess up transmit simulations
commParam.pmax       = 50;%std=40 % maximum transmission power for a bit
commParam.pmin       = 0;  % minimum transmission power for a bit
commParam.pathloss   = 1./calcPathLoss([sensorState.distToFC],2);
commParam.PtotConstr = numSensors/5*BW*commParam.pmax/log(1/commParam.PrBitError); % or numSensors*dim*2*pmax*0.8;
commParam.Wp         = Wpos; % dynamic range of position of source (*to be adjusted*)
commParam.Wv         = Wvel; % dynamic range of velocity of source (*to be adjusted*)
commParam.DynRgStrPos= [scenario.observRegion(1); scenario.observRegion(3)]; % the dynamic range starting position
commParam.fd         = dopplerFreqMax;
commParam.h          = 1*ones(numSensors,Horizon); % fading coefficients
commParam.hactual    = 1*ones(numSensors,Horizon); % fading coefficients - for when the optimization values are based on predictions 
commParam.c          = 20; % receiver and coding gain
commParam.hHistoryOracle = []; % fading coefficient history
commParam.hHistory.prevH = []; % helper variables
commParam.hHistory.time = [];
commParam.W = repmat(...
           [ones(dim,1)*commParam.Wp;...
            ones(dim,1)*commParam.Wv],...
            numSensors,1); % make a dynamic range vector that is appropriate for pos/vel different ranges
commParam.WoneSensor = ...
           [ones(dim,1)*commParam.Wp;...
            ones(dim,1)*commParam.Wv];
% Store the comm parameters
scenario.commParam = commParam;

%% Set up the energy harvester parameters
energyBoundFunc = @(W,P,dNm1,dN,T,K) W^2 / P .* (1./sqrt(dNm1) - 1./sqrt(dN))./...
                (1./(dNm1.*sqrt(dN)) - 1./(dN.*sqrt(dNm1))) * T / 15;% K = 15 : using largest K amongst variation to show improvement for increasing no. of time steps;
y = sort([commParam.pathloss 1])*log(1./commParam.PrBitError)/commParam.c;
maxAndNextMaxLoss = y(end-1:end);

scenario.energyHarvester.dynamicBatteryDet = 0;
if scenario.energyHarvester.dynamicBatteryDet
    batteryMax   = energyBoundFunc(max(commParam.W),mean([scenario.p0 scenario.q]),maxAndNextMaxLoss(2),maxAndNextMaxLoss(1),scenario.tsample,nsamples);
else
    switch scenario.objectiveType
    case 'BLUEwfwf'
        batteryMax   = 5;%energyBoundFunc(max(commParam.W),mean([scenario.p0 scenario.q]),maxAndNextMaxLoss(2),maxAndNextMaxLoss(1),scenario.tsample,nsamples);
    case 'sumrate'
        batteryMax   = 50;%energyBoundFunc(max(commParam.W),mean([scenario.p0 scenario.q]),maxAndNextMaxLoss(2),maxAndNextMaxLoss(1),scenario.tsample,nsamples);
    end
end

switch scenario.objectiveType
    case 'BLUEwfwf'
        % Temporarily using inf battery to test viable of method for batteryless system
        scenario.energyHarvester.batteryMax   = inf; % batteryMax;
        
        scenario.energyHarvester.maxMinEnergy = [.1 .4]*batteryMax; % in eu units
        scenario.energyHarvester.arrivalRateList = [1 1.5 2]'; % number of energy arrivals per time period tsamplei
    case 'sumrate'
        scenario.energyHarvester.batteryMax   = batteryMax; % batteryMax;
        scenario.energyHarvester.maxMinEnergy = [.04 .1]*batteryMax;%[.1 .3]*batteryMax; % in eu units
        scenario.energyHarvester.arrivalRateList = [3 7 10]'; % number of energy arrivals per time period tsample
    otherwise
        error('SimParams:BatterySetup','need to define a battery setup for this problem type')
end
scenario.energyHarvester.arrivalRateMCtransition = eye(3);...[.899 .1 .001;...
                                                   ... .1   .8 .1;...
                                                   ... .001  .2 .799];


rateSelect = repmat([1 2 3],1,ceil(numSensors/3));
rateSelect = rateSelect(1:numSensors);
% rateSelect = ones(1,numSensors);
%rateSelect = ceil(length(scenario.energyHarvester.arrivalRateList)*rand(1,numSensors));
arrivalRates = [];
for ii = 1:numSensors
    arrivalRates(ii) = scenario.energyHarvester.arrivalRateList(rateSelect(ii));%#ok
end
scenario.energyHarvester.arrivalRate = arrivalRates; % arrivals/second -> vector of current arrival rates

% Set the initial battery level
scenario.PowerInit  = .3*batteryMax;
scenario.PwrInitIndiv = nodeEnergyInitScale*scenario.PowerInit; % the power remaining for each of the sensor nodes


%% Set up the global sensorData variable
sensorData.commParam  = commParam;
sensorData.numSensors = numSensors;
sensorData.energyRem  = [];% will be set as <- nodeEnergyInitScale*scenario.PowerInit; % the power remaining for each of the sensor nodes
sensorData.batteryStateOptmz = zeros(numSensors,Horizon);
sensorData.numElem   = dim*2;
sensorData.useDimen    = scenario.useDimen;
sensorData.pMetricSize = []; % set later on \/ helps to determine the size of the state vector or metric being used for sensor estimate quality 
sensorData.horizon   = Horizon;
sensorData.Emax      = scenario.energyHarvester.batteryMax;
sensorData.actEnergyHarvest = []; % this vector tracks energy up current time
sensorData.energyHarvestForOptimz = [sensorData.energyRem zeros(numSensors,Horizon-1)];% current energy level and prediction of future harvest through H+k time step
sensorData.eHistoryOracle = []; % this is the omniscient pre-determined energy values (not necessarily known by the system
sensorData.pHistoryOracle = []; % state covariance history
sensorData.PstateForOptimz = [];
sensorData.PstateObjCalc   = [];
sensorData.weights = ones(numSensors,scenario.nsamples);%repmat(1:scenario.nsamples,numSensors,1); % create a future time weighted priority for log sum utilities
sensorData.useOracle = scenario.oracleMode;
sensorData.plotOptmzProgress = scenario.plotOptmzProgress;
sensorData.currentInterval = [];
sensorData.numEpochs = [];
sensorData.epochEndIdx = [];

% sensorData.objType = scenario.objectiveType;
% sensorData.Lambdas   = ones(numSensors,1);
% sensorData.mseLevel  = scenario.mseLevel;
% sensorData.mseIter   = []; % save sub-variable for use of Taylor approx of mse term <- *might need to modify this starting location with a more apt choice*


%% ============== Optimization Setup and  Variable Management =============
% Form and solve problem to optimally estimate the state from the state
% vectors reported by the sensors. First specify some options.
scenario.options = optimset('fmincon');
% options.LargeScale = 'off';
scenario.options.Algorithm = 'interior-point';%'active-set';%
scenario.options.Display = 'notify';%'iter'

%options = optimset('Display','notify'); %'iter'
scenario.options.TolFun      = 1e-6;
scenario.options.TolCon      = 1e-4;
scenario.options.TolX        = 1e-6;
scenario.options.MaxIter     = 3000;
scenario.options.MaxFunEvals = 1e6;

if scenario.singleShot
    scenario.options.MaxIter     = 3500;
    scenario.options.MaxFunEvals = 5e5;
end

% I think these two control options below are determined by the
% size of the problem, and since that's what I wanted to put down,
% I'll just leave them for now.
%   options.MaxFunEvals = 2000;
%   options.MaxIter = 2000;

scenario.options.Diagnostics = 'off';
%options.PlotFcns = {@optimplotfval,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderop});

% Options for adjusting the Sequential Convex Programming (SCP) method and 
% feasible initial vector searching
scenario.scpOptions.objChangeTol        = 1e-4;
scenario.scpOptions.maxIter             = 1; % linear upper is a single iteration solve (more iterations for methods req. mutliple SCP iterations)
scenario.scpOptions.maxIterFeasibleFind = 10;
scenario.scpOptions.powerReductionStep  = 2; % when attempting to obtain a feasible initial point, the amount by which to reduce the power 

% ==== Parameterized Tests (sensitivity to variables, unknowns, etc.) =====

scenario.tests.covSensitivity.pertbVar = 0;


%==========================================================================

% Define which functions should run for simulation and optimization
scenario.optimizerFnc              = optimizerFnc;

scenario.optFormSetup.horizon      = Horizon;
scenario.optFormSetup.simulate     = @simulate_eh;
scenario.optFormSetup.postOptmCalc = @calcPostOptmQuantitiesLinUB;

if scenario.oracleMode
    scenario.predictRoutine = @predictAllTruth;
else
    scenario.predictRoutine = @predictAllEstimate;
end

% optimization/metric functions
switch scenario.objectiveType
    case 'BLUEwfwf'
        scenario.usesObjectTraj = 1;
        scenario.optFormSetup.objFunc         = @ObjFunc_eh_wf;
        scenario.optFormSetup.constrFunc      = @constraints_eh_DE_lin;
        %scenario.optFormSetup.JcostFunc       = @JcostCompare;%@(p,q,c)1/(sum(1./(p+q+c))); % for comparison with true problem formulation
        scenario.optFormSetup.extractBitsFunc = @(dynLevels,h) log2((dynLevels+1));%@(dynLevels) log2((sqrt(dynLevels)+1));
        scenario.optFormSetup.energyUsed = @(q,sensorData)q.*reshape( repmat(...
                                                                        reshape(...
                                                                        repmat(sensorData.commParam.pathloss',1,sensorData.numEpochs)...
                                                                        ./sensorData.commParam.h...
                                                                        ,1,sensorData.numEpochs*sensorData.numSensors)...
                                                                        ,sensorData.pMetricSize,1),...
                                                                        ...
                                                                        sensorData.numSensors*sensorData.pMetricSize,sensorData.numEpochs)...
                                                               *log(1./sensorData.commParam.PrBitError)...
                                                               /sensorData.commParam.c;
        if any(strcmp(scenario.vectorMetricType,{'onedimen','trace'}))
            scenario.optFormSetup.perEpochCostFunc = @(q,sensorData)...
                                                -sum(q(:,1)./(q(:,1).*...
                                                sensorData.PstateForOptimz(:,1)+...
                                                repmat(sensorData.commParam.W(sensorData.useDimen).^2,length(q(:,1)),1)));% needs to be updated for multi-dimensional
        else
            scenario.optFormSetup.perEpochCostFunc = @(q,sensorData)...
                                                -sum(q(:,1)./(q(:,1).*...
                                                sensorData.PstateForOptimz(:,1)+...
                                                sensorData.commParam.W.^2));% needs to be updated for multi-dimensional
        end
        scenario.optFormSetup.margUtilFunc    = @(r,W) W.^2./(r+W.^2).^2;
        scenario.optFormSetup.margUtilFuncInv = @(lambda,W)(W./sqrt(lambda)-W.^2);

        switch scenario.vectorMetricType
            case 'trace'
                sensorData.commParam.Wuse = repmat((commParam.Wp+commParam.Wv)/2,numSensors,1);
                scenario.optFormSetup.weightFunction = @(b,sensorData) 1./(sensorData.PstateForOptimz(:,1)+sensorData.commParam.Wuse.^2./(2.^(b-1)));
            case 'vector'
                sensorData.commParam.Wuse  = sensorData.commParam.W;
                scenario.optFormSetup.weightFunction = @(b,sensorData) 1./(sensorData.PstateForOptimz(:,1)+sensorData.commParam.W.^2./(2.^(b-1)));        
            case 'onedimen'
                if scenario.useDimen <= scenario.numDimen;
                    sensorData.commParam.Wuse = repmat(commParam.Wp,numSensors,1);
                else
                    sensorData.commParam.Wuse = repmat(commParam.Wv,numSensors,1);
                end
            scenario.optFormSetup.weightFunction = @(b,sensorData) 1./(sensorData.PstateForOptimz(:,1)+sensorData.commParam.Wuse.^2./(2.^(b-1)));
        end

                                            
    case 'BLUErwfrwf'
        scenario.usesObjectTraj = 1;
        scenario.optFormSetup.objFunc         = @ObjFunc_eh_rwf;
        scenario.optFormSetup.constrFunc      = @constraints_eh_simple;
        %scenario.optFormSetup.JcostFunc       = @JcostCompare;%@(p,q,c)1/(sum(1./(p+q+c))); % for comparison with true problem formulation
        scenario.optFormSetup.extractBitsFunc = @(dynLevels,h) log2((dynLevels+1));%@(dynLevels) log2((sqrt(dynLevels)+1));
        scenario.optFormSetup.energyUsed = @(q,sensorData)q.*reshape( repmat(...
                                                                          reshape(...
                                                                              repmat(sensorData.commParam.pathloss',1,sensorData.numEpochs)...
                                                                              ./sensorData.commParam.h...
                                                                              ,1,sensorData.numEpochs*sensorData.numSensors)...
                                                                          ,sensorData.pMetricSize,1),...
                                                                           ...
                                                               sensorData.numSensors*sensorData.pMetricSize,sensorData.numEpochs)...
                                                           *log(1./sensorData.commParam.PrBitError)...
                                                           ...
                                                               /sensorData.commParam.c;%.* repmat(sensorData.currentInterval,sensorData.numSensors*sensorData.pMetricSize,sensorData.numEpochs);
%  --maybe not correct?--       scenario.optFormSetup.margUtilFunc    = eval(sprintf('@(r)(r+%g)./r',commParam.Wp^2)); % needs to be updated for multi-dimensional
%                               scenario.optFormSetup.margUtilFuncInv = eval(sprintf('@(lambda)(sqrt(%g./lambda)-%g)',commParam.Wp^2,commParam.Wp^2));
        scenario.optFormSetup.perEpochCostFunc = @(q,sensorData)...
                                                -sum(q(:,1)./(q(:,1).*...
                                                sensorData.PstateForOptimz(:,1)+...
                                                repmat(sensorData.commParam.Wp,length(q(:,1)),1)));% needs to be updated for multi-dimensional
    case 'sumrate'
        scenario.usesObjectTraj = 0;
        scenario.optFormSetup.objFunc         = @ObjFunc_eh_log;
        scenario.optFormSetup.constrFunc      = @constraints_eh_lin;
        %scenario.optFormSetup.JcostFunc       = @(p,sensorData)sum(log(1+p.*sensorData.commParam.h)) ;
        scenario.optFormSetup.extractBitsFunc = @(p,h)log2(1+p.*h);
        scenario.optFormSetup.energyUsed      = @(p,sensorData) p.*repmat(sensorData.currentInterval,sensorData.numSensors*sensorData.pMetricSize,1);
        scenario.optFormSetup.layeringSolver  = @layeredDirectionalSTWF;
        
        scenario.optFormSetup.margUtilFunc    = @(r,W)1./(1+r);
        scenario.optFormSetup.margUtilFuncInv = @(lambda,W)(1./lambda) - 1;
        scenario.optFormSetup.perEpochCostFunc = @(p,sensorData) -sum(log(1+p(:,1).*sensorData.commParam.h(:,1)));
    case 'weightedsumrate'
        scenario.optFormSetup.objFunc         = @ObjFunc_eh_wtlog;
        scenario.optFormSetup.constrFunc      = @constraints_eh_lin;
        %scenario.optFormSetup.JcostFunc       = @(p,sensorData)sum(log(1+p.*sensorData.commParam.h)) ;
        scenario.optFormSetup.extractBitsFunc = @(p,h)log2(1+p.*h);
        scenario.optFormSetup.energyUsed      = @(p,sensorData) p.*repmat(sensorData.currentInterval,sensorData.numSensors*sensorData.pMetricSize,1);
        
        scenario.optFormSetup.margUtilFunc    = @(r,W)1./(1+r);
        scenario.optFormSetup.margUtilFuncInv = @(lambda,W)(1./lambda) - 1;
        scenario.optFormSetup.perEpochCostFunc = @(p,sensorData) -sum(sensorData.weights.*log(1+p(:,1).*sensorData.commParam.h(:,1)));
end
switch scenario.vectorMetricType
    case 'trace'
        scenario.optFormSetup.PmetricFunc   = @(x)trace(x);
    case 'vector'
        scenario.optFormSetup.PmetricFunc   = @(x)diag(x)';
    case 'onedimen'
        scenario.optFormSetup.PmetricFunc   = @(x)x(scenario.useDimen,scenario.useDimen);
end


%--------------------------------------------------------------------------
% Define functions for extracting the appropriate information from the
% optimization vector - only have one formulation right now
%--------------------------------------------------------------------------
%scenario.optFormSetup.JcostFunc       = @JcostCompare;%@(p,q,c)1/(sum(1./(p+q+c))); % for comparison with true problem formulation

% scenario.optFormSetup.powerUsed = @(power,bits,numElem)sum(reshape(power.*bits,numElem,length(power)/numElem),1)';
%scenario.optFormSetup.extractBitsFunc = @(dynLevels) log2((dynLevels+1));% @(dynLevels) log2((sqrt(dynLevels)+1));
                                        

if isscalar(scenario.optFormSetup.PmetricFunc(eye(2)))
    scenario.optFormSetup.pMetricSize     = 1;
    scenario.optFormSetup.numVariables    = Horizon*numSensors;
    scenario.optFormSetup.numConstraints  = 1 + (Horizon+1)*numSensors;

    sensorData.pMetricSize = 1; % helps to determine the size of the state vector or metric being used for sensor estimate quality 
else
    scenario.optFormSetup.pMetricSize     = numElem;
    scenario.optFormSetup.numVariables    = Horizon*numSensors*numElem;
    scenario.optFormSetup.numConstraints  = 1 + (Horizon+1)*numSensors*numElem;
    
    sensorData.pMetricSize = numElem; % helps to determine the size of the state vector or metric being used for sensor estimate quality
end
scenario.optFormSetup.x0Func = @x0FuncLinUB;
   
scenario.optFormSetup.AeqFunc=@(numActive,commParam)[]; 
scenario.optFormSetup.beqFunc=@(numActive,commParam)[];
scenario.optFormSetup.AFunc=@(numActive,commParam)[];  
scenario.optFormSetup.bFunc=@(numActive,commParam)[];
scenario.optFormSetup.lbFunc=@(numActive,commParam)[]; 
scenario.optFormSetup.ubFunc=@(numActive,commParam)[];
