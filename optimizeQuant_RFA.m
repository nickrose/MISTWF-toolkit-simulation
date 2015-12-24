function [bits,energyUsed, Jcost, JcostThisEpoch, exitflag] = optimizeQuant_RFA(sensorData,scenario)
% optimizeQuant_RFA.m - a comparison method, randomly assigns resources
% while maintaining feasibility.
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

global recycleMCdraws

% -------------------------------------------- \/
activeSensors = ... % all sensors (excludes inoperable): 1 3 ... 7 ... N
    find(sensorData.energyRem > sensorData.commParam.pmin); % updates which sensors are active in absolute reference
numActive = length(activeSensors); % number active

q0 = zeros(numActive*scenario.optFormSetup.pMetricSize*sensorData.horizon,1);
    %= 40*ones(numActive*scenario.optFormSetup.pMetricSize*sensorData.horizon,1);%findInitialVector(tNowIdx,activeNow,numElem,optFormSetup,...numSensors,optFormSetup,...
%  reportsAtFC,scpOptions);%,scenario.reuseScenarioSetup,scenario.runIdx);
tSteplength = repmat(sensorData.currentInterval,numActive*scenario.optFormSetup.pMetricSize,1);
tSteplength = tSteplength(:);

q = q0;
Aeq = [];   beq = [];
A   = [];   b   = [];
lb  = [];   ub  = [];

vecLength = length(q0);
untriedSet = 1:vecLength;
tryNo = 1;
rvSample = fixedSampleDatabase('RFAalgSelect',@rand,recycleMCdraws.reuseNoise,tryNo,recycleMCdraws.runIdx);
[a,tryIdx] = min(abs(round(vecLength*rvSample) - untriedSet));
% iter = 0; % for debug
while ~isempty(untriedSet)
    q(tryIdx) = q(tryIdx) + 1e3;
    violMax = max(scenario.optFormSetup.constrFunc(q));
    if violMax > 0
        q(tryIdx) = q(tryIdx) - violMax/tSteplength(tryIdx);
        untriedSet = setdiff(untriedSet,tryIdx);

        tryNo = tryNo + 1;
        rvSample = fixedSampleDatabase('RFAalgSelect',@rand,recycleMCdraws.reuseNoise,tryNo,recycleMCdraws.runIdx);
        [a,idx] = min(abs(round(length(untriedSet)*rvSample) - untriedSet));
        tryIdx = untriedSet(idx);
    end
%     iter = iter + 1;
%     if iter > 10000 && iter < 10002
%         keyboard;
%     end
end
q(q<0) = 0;

Jcost = scenario.optFormSetup.objFunc(q);
exitflag = 1;
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

% TEMP: just need for globecom paper, will figure out better way to
% make this general later
JcostThisEpoch = -sum(log(1+perSensorResourceLevels(:,1).*sensorData.commParam.h(:,1)));





