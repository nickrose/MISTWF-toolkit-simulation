function [bits,energyUsed, Jcost, JcostThisEpoch, exitflag] = optimizeQuant_layered(sensorDataIn,scenario)
% optimizeQuant_layered.m - a simple way to contain the homebrew 'layered'
% solution, simple instaintiation and control of variable names, etc.
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
% By: Nick Roseveare, Feb 2012
global sensorData

% -------------------------------------------- \/
activeSensors = ... % all sensors (excludes inoperable): 1 3 ... 7 ... N
    find(sensorData.energyRem >= 0);%sensorData.commParam.pmin); % updates which sensors are active in absolute reference
numActive = length(activeSensors); % number active

misc.plotProgress = sensorData.plotOptmzProgress;
misc.playBack = 0;
switch scenario.objectiveType
    case 'BLUEwfwf'
        params.a = 1./sensorData.PstateForOptimz;
        params.b = sensorData.PstateForOptimz;
        params.c = zeros(numActive*sensorData.pMetricSize,sensorData.numEpochs);
        d = (reshape( repmat(...
                                reshape(...
                                    repmat(sensorData.commParam.pathloss',1,sensorData.numEpochs)...
                                    ./sensorData.commParam.h...
                                ,1,sensorData.numEpochs*numActive)...
                            ,sensorData.pMetricSize,1),...
                        ...
                   numActive*sensorData.pMetricSize,sensorData.numEpochs)...
                    *log(1./sensorData.commParam.PrBitError)...
                    /sensorData.commParam.c...
                    ).^2; % <- the alpha*h terms are squared since the quantz level variable is a squared result
        
        % the scaling cost of resources for actual fading and state estimate variance
        dactual = (reshape( repmat(...
                                reshape(...
                                    repmat(sensorData.commParam.pathloss',1,sensorData.numEpochs)...
                                    ./sensorData.commParam.hactual...
                                ,1,sensorData.numEpochs*numActive)...
                            ,sensorData.pMetricSize,1),...
                        ...
                   numActive*sensorData.pMetricSize,sensorData.numEpochs)...
                    *log(1./sensorData.commParam.PrBitError)...
                    /sensorData.commParam.c...
                    ).^2;
                
        % including redefinition of r, b, d for modified form:
        params.b = params.b ./ d;
        params.d = ones(size(d));
        
        misc.W = sensorData.commParam.Wuse; 
        E = (sensorData.energyHarvestForOptimz./repmat(sensorData.currentInterval,sensorData.numSensors,1)).^2;
        PtotMax = sensorData.commParam.PtotConstr^2;
        Emax = sensorData.Emax^2;
        
    case 'sumrate'
        params.a = ones(numActive,sensorData.numEpochs);
        params.b = 1./sensorData.commParam.h;
        params.c = zeros(numActive,sensorData.numEpochs);
        params.d = ones(numActive,sensorData.numEpochs);
        
        params.b = params.b ./ params.d;
        params.d = ones(size(params.d));
        
        E = sensorData.energyHarvestForOptimz;
        PtotMax = sensorData.commParam.PtotConstr;
        Emax = sensorData.Emax;
    case 'weightedsumrate'
        params.a = sensorData.weights;
        params.b = 1./sensorData.commParam.h;
        params.c = zeros(numActive,sensorData.numEpochs);
        params.d = ones(numActive,sensorData.numEpochs);
        
        params.b = params.b ./ params.d;
        params.d = ones(size(params.d));
        
        E = sensorData.energyHarvestForOptimz;
        PtotMax = sensorData.commParam.PtotConstr;
        Emax = sensorData.Emax;
end

[r,exitflag] = generalMISTWF(scenario.optFormSetup.margUtilFunc,scenario.optFormSetup.margUtilFuncInv,params,...
                                E,Emax,PtotMax,misc);
                            
                            % including redefinition of r, b for modified
                            % form
sensorData.PstateForOptimz = sensorData.PstateObjCalc; % change the estimate uncertainty to the actual uncertainty vs. the predicted uncertainty
if strcmp(scenario.objectiveType,'BLUEwfwf')
    Jcost = scenario.optFormSetup.objFunc(r(:)./dactual(:));
else
    Jcost = scenario.optFormSetup.objFunc(r(:));
end

if exitflag < 0
    disp('no feasible point found')
    %keyboard
elseif exitflag == 0
    disp('feasible: not enough iterations')
else
    disp('solution found')
end

% -------------------------------------------- /\
numEpochs = sensorData.numEpochs;

if strcmp(scenario.objectiveType,'BLUEwfwf')
    % including redefinition of r, b, (by d) for modified form
    perSensorResourceLevelsTrueForObj = reshape(r./dactual,scenario.optFormSetup.pMetricSize*numActive,numEpochs);
    perSensorResourceLevels = reshape(r./d,scenario.optFormSetup.pMetricSize*numActive,numEpochs);
else
    perSensorResourceLevels = reshape(r,scenario.optFormSetup.pMetricSize*numActive,numEpochs);
end
bits = scenario.optFormSetup.extractBitsFunc(perSensorResourceLevels,sensorData.commParam.h);

energyUsed = scenario.optFormSetup.energyUsed(perSensorResourceLevels,sensorData);

sensorData.batteryStateOptmz = zeros(numActive,numEpochs);
% Battery not updated in the layered wf function; must do it here
for jj = 1:numEpochs
    sensorData.batteryStateOptmz(:,jj) = sum(sensorData.energyHarvestForOptimz(:,1:jj) - energyUsed(:,1:jj),2);
end

% TEMP: just need for globecom paper, will figure out better way to
% make this general later
if any(strcmp(scenario.vectorMetricType,{'onedimen','trace'}))
    if strcmp(scenario.objectiveType,'BLUEwfwf')
        JcostThisEpoch = scenario.optFormSetup.perEpochCostFunc(perSensorResourceLevelsTrueForObj,sensorData);
    else
        JcostThisEpoch = scenario.optFormSetup.perEpochCostFunc(perSensorResourceLevels,sensorData);
    end
else
    if strcmp(scenario.objectiveType,'BLUEwfwf')
        JcostThisEpoch = scenario.optFormSetup.perEpochCostFunc(r./dactual,sensorData);
    else
        JcostThisEpoch = scenario.optFormSetup.perEpochCostFunc(r,sensorData);
    end
end

