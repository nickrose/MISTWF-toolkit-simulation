% function runEHscenarios()
% runEHscenarios.m - Code for running multiple Energy Harvesting WSN tests.
%
% This script simulates the energy harvesting WSN with the application of
% distributed estimation or sum log rate.
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
% By: Nick Roseveare, Dec 2011
%clear all global

global recycleMCdraws %progress

% warning('error','MATLAB:divideByZero');

%logFileID = fopen('./logs_EHscenarios.log','w');
try
    !rm logs_EHscenarios.log
end
diary('logs_EHscenarios.log')

dateStart = datestr(now);
fprintf('\n\n\n++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\nStarting date/time: %s\n',dateStart)
fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++\n\n\n')


% Variables to vary over:
varyOverList = {'nsamples','numSensors','dopplerFreqMax'};% % can be 'numSensors' or 'nsamples'
nsamplesList = 20:-round((20-6)/7):6;
numSensorsList = 20:-round((20-6)/7):6;
dopplerFreqMaxList = 2.1:-(2.1-.1)/7:.1;


numSenDefault = 5;
nSamplesDefault = 5;
dopplerFreqMaxDefault = 1;

numMC = 1;

recycleMCdraws.reuseNoise = 1;
recycleMCdraws.runIdx    = 1;

inputParams.numSensors   = 5;
inputParams.nsamples     = 10;
    Horizon = 5;
inputParams.Horizon      = Horizon;
inputParams.BW           = 1000;
inputParams.dopplerFreqMax = 1;
inputParams.oracleMode   = 1;
inputParams.singleShot   = 1;
inputParams.stopIferror  = 1;
inputParams.objType      = 'sumrate';
% inputParams.objType      = 'BLUEwfwf';

% delete the 'saved trajectory' to start new on the simulations
% delete('lastTrueFile.mat')

prevTimeAnlys = 0;
tStart = tic;
if numMC == 1
    
    recycleMCdraws.runIdx = 1;
    
    fixedSampleDatabase('cleardb');
    % SINGLE MC RUNS
    fprintf('\nRunning MIST-WF with causal MPC\n\n')
    % layered WF optimization routine - causal
    k = 1;
    inputParams.singleShot   = 0;
    inputParams.oracleMode   = 0;
    scenario = simulationParams_eh(inputParams);
    outputDataStructBlank = makeOutputDataStructEH(scenario.numSensors,scenario.numDimen,...
        scenario.time,scenario.optFormSetup.numConstraints,scenario.optFormSetup.pMetricSize);

    scenario.optimizerFnc = @optimizeQuant_layered;
    [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
    Metrics = outputData;
    Metrics.sensorData = scenario.sensorData;
    Metrics.algName = 'layered_WF_MPC';
    %     outputData.sensorData = scenario.sensorData;
    %     Metrics.predictedLayeredWF = outputData;

    
%     fprintf('\nRunning MIST-WF with non-causal single shot\n\n')
%     % layered WF optimization routine - non causal with single shot
%     k = 3;
%     inputParams.singleShot   = 1;
%     inputParams.oracleMode   = 1;
%     scenario = simulationParams_eh(inputParams);
%     outputDataStructBlank = makeOutputDataStructEH(scenario.numSensors,scenario.numDimen,...
%                 scenario.time,scenario.optFormSetup.numConstraints,scenario.optFormSetup.pMetricSize);
% 
%     scenario.optimizerFnc = @optimizeQuant_layered;
%     [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
%     MetricsT = outputData;
%     MetricsT.sensorData = scenario.sensorData;
%     MetricsT.algName = 'layered_WF_NC';
%     Metrics(k) = MetricsT;
    %     outputData.sensorData = scenario.sensorData;
    %     Metrics.noncausalknownLayeredWF = outputData;
    
    
	fprintf('\nRunning std fmincon mtd with non-causal information\n\n')
    % std optimization routine
    k = 2;
    inputParams.singleShot   = 1;
    inputParams.oracleMode   = 1;
    scenario = simulationParams_eh(inputParams);
    outputDataStructBlank = makeOutputDataStructEH(scenario.numSensors,scenario.numDimen,...
        scenario.time,scenario.optFormSetup.numConstraints,scenario.optFormSetup.pMetricSize);
    
    scenario.optimizerFnc = @optimizeQuant_std;
    [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
    MetricsT = outputData;
    MetricsT.sensorData = scenario.sensorData;
    MetricsT.algName = 'std_optmz';
    Metrics(k) = MetricsT;
    
    diary('off')
    
    % Find single MC metric
       metrics_eh('energyHarvesting',scenario,Metrics,date,1);
    % Find comparitive MC metrics
        % metricsMC_eh(['tspEHresults_' inputParams.objType],scenario,Metrics);
else
    
    for ss = 1:length(varyOverList)
        tic;
        varyOver = varyOverList{ss};
        switch varyOver % set the 'other' scenario set variable to a default
            case 'nsamples'
                inputParams.numSensors = numSenDefault;
                inputParams.dopplerFreqMax = dopplerFreqMaxDefault;
            case 'numSensors'
                inputParams.nsamples = nSamplesDefault;
                inputParams.dopplerFreqMax = dopplerFreqMaxDefault;
            case 'dopplerFreqMax'
                inputParams.numSensors = numSenDefault;
                inputParams.nsamples = nSamplesDefault;
        end
        Metrics = [];
        Metrics.numMC = numMC;
        varyOverListMat = eval([varyOver 'List']);
        varyOverListLen = length(varyOverListMat);
        
        Metrics.varyOver = varyOverListMat;
        Metrics.varyOverStr = varyOver;
        
        for vv = 1:varyOverListLen
            
            fixedSampleDatabase('cleardb');
            
            inputParams.(varyOver) = varyOverListMat(vv);
            %inputParams.nodeEnergyInitScale = ones(inputParams.numSensors,1);  % should be values 0->1
            
            scenario = simulationParams_eh(inputParams);
            
            outputDataStructBlank = makeOutputDataStructEH(scenario.numSensors,scenario.numDimen,...
                scenario.time,scenario.optFormSetup.numConstraints,scenario.optFormSetup.pMetricSize);
            
            
            
           
            tryAgain = 1;
            startRun = 1;
            while tryAgain
                try
                    
                    % ******* Method: MISTWF - MPC *******
                    
                    inputParams.singleShot   = 0;
                    inputParams.oracleMode   = 0;
                    inputParams.Horizon      = Horizon;%scenario.nsamples;
                    for mm = startRun:Metrics.numMC;
                        recycleMCdraws.runIdx = mm;
                        clc
                        fprintf('\nPrev algorithm solution took %g minutes\n\n',prevTimeAnlys);
                        
                        fprintf('\nUsing MISTWF method with causal knowledge for\n') %<<<<<<<<<<<<<<
                        fprintf('Scenario varying over %s, %d of %d scenarios\n',varyOver,vv,varyOverListLen)
                        fprintf('Running MC run %d of %d\n',mm,Metrics.numMC)
                        if scenario.energyHarvester.dynamicBatteryDet
                            fprintf('Maximum battery dynamically determined to be: %g\n',scenario.energyHarvester.batteryMax)
                        else
                            fprintf('Maximum battery set to: %g\n',scenario.energyHarvester.batteryMax)
                        end
                        % MISTWF optimization routine - causal stochastic knowledge
                        scenario = simulationParams_eh(inputParams);
                        scenario.optimizerFnc = @optimizeQuant_layered;
                        
                        
                        fprintf('\nPreliminaries took %g minutes\n',toc(tStart)/60);
                        
                        
                        %                 tryAgain = 1;
                        %                 while tryAgain % should only need this 'try-catch' on the first run when trajectories are generated for each MC run...afterward the track should be filter fine for the rest of the methods
                        %                     try
                        tStart = tic;
                        [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
                        prevTimeAnlys = toc(tStart)/60;
                        %                         tryAgain = 0;
                        %                     catch messageID
                        %                         if ~(strcmpi(messageID.identifier,'MATLAB:posdef') && length(messageID.stack) == 5)
                        %                             rethrow(messageID);
                        %                         else
                        %                             inputParams.newWSNposEachMCrun = 1;
                        %                             inputParams.newTrajEachMCrun   = 1;
                        %                             scenario = simulationParams_eh(inputParams);
                        %                         end
                        %                     end
                        %                 end
                        %                 inputParams.newWSNposEachMCrun = 0;
                        %                 inputParams.newTrajEachMCrun   = 0;
                        
                        tStart = tic;
                        outputData.sensorData = scenario.sensorData;
                        Metrics.predictedLayeredWF(mm,vv) = outputData;
                    end
                    
                    % ******* Method: MISTWF - NC *******
                    inputParams.singleShot   = 1;
                    inputParams.oracleMode   = 1;
                    for mm = startRun:Metrics.numMC;
                        recycleMCdraws.runIdx = mm;
                        clc
                        fprintf('\nPrev algorithm solution took %g minutes\n\n',prevTimeAnlys);
                        
                        fprintf('\nUsing MISTWF method with non-causal knowledge for\n')%<<<<<<<<<<<<<<
                        fprintf('Scenario varying over %s, %d of %d scenarios\n',varyOver,vv,varyOverListLen)
                        fprintf('Running MC run %d of %d\n',mm,Metrics.numMC)
                        if scenario.energyHarvester.dynamicBatteryDet
                            fprintf('Maximum battery dynamically determined to be: %g\n',scenario.energyHarvester.batteryMax)
                        else
                            fprintf('Maximum battery set to: %g\n',scenario.energyHarvester.batteryMax)
                        end
                        % MISTWF optimization routine - non-causal stochastic knowledge
                        scenario = simulationParams_eh(inputParams);
                        
                        scenario.optimizerFnc = @optimizeQuant_layered;
                        
                        
                        fprintf('\nPreliminaries took %g minutes\n',toc(tStart)/60);
                        
                        tStart = tic;
                        [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
                        prevTimeAnlys = toc(tStart)/60;
                        
                        tStart = tic;
                        outputData.sensorData = scenario.sensorData;
                        Metrics.noncausalknownLayeredWF(mm,vv) = outputData;
                    end
                    
                    % ******* Method - RFA *******
                    inputParams.singleShot   = 1;
                    inputParams.oracleMode   = 1;
                    for mm = startRun:Metrics.numMC
                        recycleMCdraws.runIdx = mm;
                        clc
                        fprintf('\nPrev algorithm solution took %g minutes\n\n',prevTimeAnlys);
                        
                        fprintf('\nUsing RFA method for\n')%<<<<<<<<<<<<<<
                        fprintf('Scenario varying over %s, %d of %d scenarios\n',varyOver,vv,varyOverListLen)
                        fprintf('Running MC run %d of %d\n',mm,Metrics.numMC)
                        if scenario.energyHarvester.dynamicBatteryDet
                            fprintf('Maximum battery dynamically determined to be: %g\n',scenario.energyHarvester.batteryMax)
                        else
                            fprintf('Maximum battery set to: %g\n',scenario.energyHarvester.batteryMax)
                        end
                        % RFA optimization routine
                        scenario = simulationParams_eh(inputParams);
                        scenario.optimizerFnc = @optimizeQuant_RFA;
                        
                        
                        fprintf('\nPreliminaries took %g minutes\n',toc(tStart)/60);
                        
                        tStart = tic;
                        [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
                        prevTimeAnlys = toc(tStart)/60;
                        
                        tStart = tic;
                        outputData.sensorData = scenario.sensorData;
                        Metrics.rfaRoutine(mm,vv) = outputData;
                    end
                    
                    % ******* Method - fmincon (NC)*******
                    inputParams.singleShot   = 1;
                    inputParams.oracleMode   = 1;
                    for mm = startRun:Metrics.numMC
                        recycleMCdraws.runIdx = mm;
                        clc
                        fprintf('\nPrev algorithm solution took %g minutes\n\n',prevTimeAnlys);
                        
                        fprintf('\nUsing std fmincon method (non-causal) for\n')%<<<<<<<<<<<<<<
                        fprintf('Scenario varying over %s, %d of %d scenarios\n',varyOver,vv,varyOverListLen)
                        fprintf('Running MC run %d of %d\n',mm,Metrics.numMC)
                        if scenario.energyHarvester.dynamicBatteryDet
                            fprintf('Maximum battery dynamically determined to be: %g\n',scenario.energyHarvester.batteryMax)
                        else
                            fprintf('Maximum battery set to: %g\n',scenario.energyHarvester.batteryMax)
                        end
                        % std optimization routine - noncausal knowledge of statistics
                        scenario = simulationParams_eh(inputParams);
                        scenario.optimizerFnc = @optimizeQuant_std;
                        
                        
                        fprintf('\nPreliminaries took %g minutes\n',toc(tStart)/60);
                        
                        tStart = tic;
                        [outputData, scenario] = simulate_eh(scenario, scenario.sensorState, outputDataStructBlank);
                        prevTimeAnlys = toc(tStart)/60;
                        
                        tStart = tic;
                        outputData.sensorData = scenario.sensorData;
                        Metrics.stdOptimizRoutine(mm,vv) = outputData;
                    end
                    
                    Metrics.scenario(vv) = scenario;
                    tryAgain = 0;
                catch messageID
                        if   ...strcmpi(messageID.identifier,'MATLAB:divideByZero')&&(length(messageID.stack) == 5)||...
                            ~(strcmpi(messageID.identifier,'MATLAB:posdef')) ||...
                            ~((length(messageID.stack) == 5)||(length(messageID.stack) == 3))
                              %||(strcmpi(messageID.identifier,'MISTWF:badprofile') )...
                            
                            rethrow(messageID);
                        else
                            if strcmpi(messageID.identifier,'MATLAB:posdef')
                                inputParams.newWSNposEachMCrun = 1;
                                inputParams.newTrajEachMCrun   = 1;
                            end
                            scenario = simulationParams_eh(inputParams);
                            scenario.optimizerFnc = @optimizeQuant_layered;
                            
                            inputParams.newWSNposEachMCrun = 0;
                            inputParams.newTrajEachMCrun   = 0;
                            
                            startRun = recycleMCdraws.runIdx;
%                             if strcmpi(messageID.identifier,'MISTWF:badprofile')
%                                 recycleMCdraws.reuseNoise = 0;
%                                 fprintf('\nBad energy usage profile: attempting to regenerate\n\n')
%                             end
                        end
                end
            end
        end
%         diary('off')
        % Find comparitive MC metrics
        metricsMC_eh(['transCommEHresults_' inputParams.objType],scenario,Metrics);
        
    end
end
fprintf('\n\n\n++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\nStarting date/time: %s\n',dateStart)
fprintf('Finishing date/time: %s\n',datestr(now))
fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++')
fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++\n\n\n')

diary('off')
