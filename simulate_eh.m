function [outputData, scenario] = simulate_eh(varargin)
% simulate_eh.m
%
%  USAGE:
%    [outputData, scenario] = simulate_eh(scenario,sensorState,reportsAtFC)
%
% This function performs the simulation of the WSN reporting and fusing
% data.
%
%
% Copyright (C) 2012  Nick Roseveare
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
% By: Nick Roseveare, 11-2011
% Updated: 03-2012

global sensorData progress

if nargin ~=0
    scenario = varargin{1};
    sensorState = varargin{2};
    outputData = varargin{3};
    %else
    %     scenario   = simulationParams_eh;
    %     outputData = makeReportStructEH(scenario.numSensors,scenario.numDimen,...
    %         scenario.time,scenario.optFormSetup.numConstraints,scenario.optFormSetup.pMetricSize);
end
%=================================
% Assign local names of variables
% options    = scenario.options;
% scpOptions = scenario.scpOptions;

time       = 0:scenario.tsample:(scenario.tsample*(scenario.nsamples-1));
numElem    = scenario.numDimen*2;
numSensors = scenario.numSensors;
numElem    = scenario.numElem;
% dfltConstr = inf*ones(scenario.optFormSetup.numConstraints,1);
optFormSetup = scenario.optFormSetup;

% Ensure that the original battery state is restored in case running for
% multiple algorithms without rerunning the parameters file
sensorData.energyRem  = scenario.PwrInitIndiv;

% disp('Before generation')
% sensorData.commParam.h

%==========================================================================
% Generate object movement data and measurements for this sim
if scenario.usesObjectTraj
    sensorStateOracle = sensorState;
    for tt = 1:scenario.nsamples
        for nn = 1:numSensors
            [z, R, t] = generateMeasurement_rangeAng(sensorStateOracle(nn),scenario.trueSys,tt,scenario.useRealisticRngAng);
            sensorStateOracle(nn) = filterMeasurement(sensorStateOracle(nn),z,R,t);
            % Save the Oracle data
            if isscalar(optFormSetup.PmetricFunc(eye(2)))
                sensorData.pHistoryOracle(nn,tt) = optFormSetup.PmetricFunc(sensorStateOracle(nn).P); % state covariance history
            else
                sensorData.pHistoryOracle(1+numElem*(nn-1):(numElem*nn),tt) = optFormSetup.PmetricFunc(sensorStateOracle(nn).P); % state covariance history
            end
            outputData.x(:,nn,tt)       = sensorStateOracle(nn).x;
            outputData.Pmetric(nn,:,tt)      = optFormSetup.PmetricFunc(sensorStateOracle(nn).PHist(:,:,tt));
            outputData.indivSensorMSE(nn,tt) = norm(scenario.trueSys.xtrue(:,tt)-sensorStateOracle(nn).xHist(:,tt));
            outputData.Pdiag(:,nn,tt)        = diag(sensorStateOracle(nn).PHist(:,:,tt)); % nonfunctioning sensors will report infinity on diagonal
        end
        outputData.P(:,:,tt) = [sensorStateOracle.P]; % take covariance as is for now (*can add noise to covariance report later*)
    end
end


% Generate fading
[R,corrSign] = retrieveGaussianCorrFromRayl(time,time,scenario.commParam.fd);

idxStr.setLen = scenario.nsamples;
for nn = 1:numSensors
    idxStr.curIdx = nn;
    h = generateCorrRayleighRVs(R,corrSign,idxStr);
    sensorData.commParam.hHistoryOracle(nn,:) = h'; % fading coefficient history
end
% disp('After generation')
% sensorData.commParam.h

% Generate energy arrival
sensorData = generateHarvestedEnergy(scenario,sensorData,time);

%==========================================================================

if scenario.singleShot % extend the horizon to the end of the scenario run and attempt to estimate all levels at once
    disp(['Running Optimization in single shot mode'])
    
    
    % Execute predictions from current time step up to horizon
    sensorData.epochEndIdx = sensorData.horizon;
    sensorData.numEpochs   = sensorData.horizon;
    sensorData = scenario.predictRoutine(sensorData,sensorState,scenario,1,time);
    
    if scenario.plotOptmzProgress
        progress = resetProgress(numSensors,sensorData.horizon,scenario);
    end
    
    [bits,energyUsed, Jc, JcFirstEpoch, exitflag] = scenario.optimizerFnc(sensorData,scenario);
    
    bitsPerSensor = reshape(sum(reshape(bits,scenario.optFormSetup.pMetricSize,numSensors*sensorData.numEpochs),1),numSensors,sensorData.numEpochs);
    energyUsedPerSensor = reshape(sum(reshape(energyUsed,optFormSetup.pMetricSize,numSensors*sensorData.numEpochs),1),numSensors,sensorData.numEpochs);
    
    outputData.objectiveVal               = Jc;
    %outputData.OptObjCompare             = Jd;
    outputData.exitCode                   = exitflag;
    outputData.bits                       = bits;
    outputData.bitsPerSensor              = bitsPerSensor;
    outputData.energyUsed                 = energyUsed;
    outputData.energyUsedPerSensor        = energyUsedPerSensor;
    outputData.batteryState               = sensorData.batteryStateOptmz;
    
else
    %% Run time-stepping Kalman filter and optimally estimate state vector
    disp(['Running Optimization for possibly ', num2str(scenario.nsamples),' different time steps'])
    
    % reset the sensor states to beginning (since sensors could be
    % turned off, the increase of error of a local state vector could
    % increase
    fprintf('\n');
    for tNowIdx = 1:scenario.nsamples
        fprintf('Processing time index: %d ',tNowIdx)
      
        if scenario.usesObjectTraj
            % Get the newest measurement and update the state
            for nn = 1:numSensors
                [z, R, t] = generateMeasurement_rangeAng(sensorState(nn),scenario.trueSys,tNowIdx,scenario.useRealisticRngAng );
                sensorState(nn) = filterMeasurement(sensorState(nn),z,R,t);
                % Oracle data saved at beginning of simulation
            end
        end
        %==================================================================
        % Execute predictions from current time step up to horizon (or end
        % of the scenario)
        sensorData.epochEndIdx = min(tNowIdx-1+sensorData.horizon,scenario.nsamples);
        sensorData.numEpochs   = min(sensorData.horizon,scenario.nsamples-tNowIdx+1);

        sensorData = scenario.predictRoutine(sensorData,sensorState,scenario,tNowIdx,time);
        

            %         disp('After prediction')
            %         sensorData.commParam.h
        
        % Find the optimal power allocations through the horizon
        %==========\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=\/=====================
        [bits, energyUsed, Jc, JcFirstEpoch, exitflag] = scenario.optimizerFnc(sensorData,scenario);
        %==========/\=/\=/\=/\=/\=/\=/\=/\=/\=/\=/\=/\=====================
        

            %         disp('After optimization')
            %         sensorData.commParam.h
        
        bitsPerSensor = reshape(sum(reshape(bits,scenario.optFormSetup.pMetricSize,numSensors*sensorData.numEpochs),1),numSensors,sensorData.numEpochs);
        energyUsedPerSensor = reshape(sum(reshape(energyUsed,optFormSetup.pMetricSize,numSensors*sensorData.numEpochs),1),numSensors,sensorData.numEpochs);
        
        sensorData.energyRem = max(0,min(sensorData.energyRem+sensorData.eHistoryOracle(:,tNowIdx)-energyUsedPerSensor(:,1), scenario.energyHarvester.batteryMax)); % update the amt of power left at the node
        
        if scenario.usesObjectTraj
        % need for later, must finish simulation of sphere decoded QAM
        % signals
        %         thisX = transmitSeqSim([sensorState.x],activeIndex,...
        %             reportsAtFC(tNowIdx).Pmetric,bits,power,scenario.numDimen,...
        %             commParam,tNowIdx,scenario.fixNoiseSamples,scenario.runIdx);
        
        % Meanwhile, just perform the quantization of the estimate
            if strcmp(scenario.vectorMetricType,'onedimen')
                elemVec = scenario.useDimen;
                xreports = zeros(1,numSensors);
            else
                elemVec = 1:numElem;
                xreports = zeros(numElem,numSensors);
            end
            bits = max(0,min(bits,scenario.commParam.BWmax));
            for nn = 1:numSensors
                if strcmp(scenario.vectorMetricType,'onedimen')
                    if scenario.useDimen <= scenario.numDimen
                        W = sensorData.commParam.Wp;
                        regions = [(sensorData.commParam.DynRgStrPos(scenario.useDimen)+W/(2^bits(nn,1))):(2*W/(2^bits(nn,1))):2*W];% middle of regions - toward middle of each
                    else
                        W = sensorData.commParam.Wv;
                        regions = [(-W+W/(2^bits(nn,1))):(2*W/(2^bits(nn,1))):W];% middle of regions - toward middle of each
                    end
                    [val, quantDataIdx] = min(abs(regions-outputData.x(scenario.useDimen,nn,tNowIdx)));
                                                    
                    xreports(1,nn) = regions(quantDataIdx);
                else
                    for dd = elemVec
                        if dd <=scenario.numDimen
                            W = sensorData.commParam.Wp;
                            regions = [(sensorData.commParam.DynRgStrPos(dd)+W/(2^bits(dd+numElem*(nn-1),1))):(2*W/(2^bits(dd+numElem*(nn-1),1))):2*W];% middle of regions - toward middle of each
                        else
                            W = sensorData.commParam.Wv;
                            regions = [(-W+W/(2^bits(dd+numElem*(nn-1),1))):(2*W/(2^bits(dd+numElem*(nn-1),1))):W];% middle of regions - toward middle of each
                        end
                        [val, quantDataIdx] = min(abs(regions-outputData.x(dd,nn,tNowIdx)));

                        xreports(dd,nn) = regions(quantDataIdx);
                    end
                end
            end
        
            % Find the weights and estimate based on the quantization levels
            weights  = optFormSetup.weightFunction(bitsPerSensor(:,1),sensorData);
            if strcmp(scenario.vectorMetricType,'onedimen')
                estVal   = sum(reshape(xreports(:).*weights,1,numSensors),2) ./ sum(reshape(weights,1,numSensors),2);
            else
                estVal   = sum(reshape(xreports(:).*weights,numElem,numSensors),2) ./ sum(reshape(weights,numElem,numSensors),2);
            end
        end
        % ********************************************************************
        % RECORD DATA: Finish recording fusion center data
        
        %outputData.optWeight(activeVarHelperIdx,kk) = weights;
        %outputData.q(:,tNowIdx)   = ;
        
        if scenario.usesObjectTraj
%             outputData.x(:,tNowIdx)              =
            outputData.xEstHist(:,tNowIdx)       = estVal;
            if strcmp(scenario.vectorMetricType,'onedimen')
                outputData.xError(tNowIdx)     = norm(estVal-scenario.trueSys.xtrue(scenario.useDimen,tNowIdx));
            else
                outputData.xError(tNowIdx)     = norm(estVal-scenario.trueSys.xtrue(:,tNowIdx));
            end
        end
        
        outputData.objectiveVal(tNowIdx)      = JcFirstEpoch;
        outputData.exitCode(tNowIdx)          = exitflag;
        outputData.bits(:,tNowIdx)            = bits(:,1);
        outputData.bitsPerSensor(:,tNowIdx)   = bitsPerSensor(:,1);
        outputData.energyUsed(:,tNowIdx)      = energyUsed(:,1);
        outputData.energyUsedPerSensor(:,tNowIdx) = energyUsedPerSensor(:,1);
        outputData.batteryState(:,tNowIdx)    = sensorData.energyRem;
        
        %outputData.OptObjCompare(tNowIdx)     = Jd;
        %outputData.activeIndex(:,tNowIdx)     = activeIndex';
        %outputData.activeIndex(:,tNowIdx)     = activeIndex';
        
    end
end
scenario.sensorState = sensorState; % save the sensor KF processing results
scenario.sensorData  = sensorData;

%=========================================================

function progress = resetProgress(N,H,scenario)
progress.onlinePlotting = scenario.plotOptmzProgOnline;

% progress.figures.objPlotHandle = figure(11);
% progress.figures.constrPlotHandle = figure(12);

progress.numSensors = N;
% progress.horizon = H;
progress.objFunc = scenario.optFormSetup.objFunc;

progress.obj = [];
progress.quantzLevels = [];

progress.constrValue = [];
progress.constrViol = [];

ColorOrder = ...
    [0    0     1;...
    0   .5     0;...
    1    0     0;...
    .5   0     1;...
    .75  0   .75;...
    .75 .75    0;...
    .25 .25  .25];

progress.ColorOrder = [ColorOrder; ColorOrder*diag([.8 .7 .5]); ColorOrder*diag([.5  .4 .8]) ];

