function [bits,energyUsed, Jcost, JcostThisEpoch, exitflag] = optimizeQuant_std(sensorData,scenario)
% optimizeQuant_std.m - a simple way to contain the matlab 'fmincon' based
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


% -------------------------------------------- \/
activeSensors = ... % all sensors (excludes inoperable): 1 3 ... 7 ... N
    find(sensorData.energyRem > sensorData.commParam.pmin); % updates which sensors are active in absolute reference
numActive = length(activeSensors); % number active

q0 = zeros(numActive*scenario.optFormSetup.pMetricSize*sensorData.numEpochs,1);
    %= 40*ones(numActive*scenario.optFormSetup.pMetricSize*sensorData.horizon,1);%findInitialVector(tNowIdx,activeNow,numElem,optFormSetup,...numSensors,optFormSetup,...
%  reportsAtFC,scpOptions);%,scenario.reuseScenarioSetup,scenario.runIdx);
q0=min(1,max(0,randn(length(q0),1)));
q = q0;
Aeq = [];   beq = [];
A   = [];   b   = [];
lb  = [];   ub  = [];

[q,Jcost,exitflag,output,lambda] = fmincon(scenario.optFormSetup.objFunc,q,A,b,Aeq,beq,lb,ub,scenario.optFormSetup.constrFunc,scenario.options);

if exitflag < 0
    disp('no feasible point found')
    %keyboard
elseif exitflag == 0
    disp('feasible: not enough iterations')
else
    disp('solution found')
end
%         % ***** Deal with optmz when no feasible points exist *****
%         if exitflag < 0 % no feasible points found
%             Jc = inf;
%             Jd = inf;
%             %estVal = inf*ones(2*scenario.numDimen,1);
%             weights = inf*ones(2*scenario.numDimen*numSensors,1);
%         elseif scenario.scpOptions.maxIter > 0
%             lagrMult(lagrIndex) = lambda.ineqnonlin;
%         end

% -------------------------------------------- /\
if any(any(q < -1e-4))
    disp('WARNING: constraint not working, negative "q" values selected')
end
q(q < 0) = 0;
numEpochs = sensorData.numEpochs;
perSensorResourceLevels = reshape(q,scenario.optFormSetup.pMetricSize*scenario.numSensors,numEpochs);
bits = scenario.optFormSetup.extractBitsFunc(perSensorResourceLevels,sensorData.commParam.h);

energyUsed = scenario.optFormSetup.energyUsed(perSensorResourceLevels,sensorData);

% Calculate and return the cost function in this epoch
JcostThisEpoch = scenario.optFormSetup.perEpochCostFunc(perSensorResourceLevels,sensorData);




