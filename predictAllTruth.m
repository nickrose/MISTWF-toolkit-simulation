function sensorData = predictAllTruth(sensorData,sensorState,scenario,tNowIdx,time)
% predictAllTruth.m - "predicts" the fading, energy arrivals, and state of
% the moving target. Allows for oracle-like testing of the rest of the
% code.
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
% By: Nick Roseveare, December 2011

% "Predict" sensor noise variances

sensorData.PstateForOptimz = [];
if scenario.usesObjectTraj
    sensorData.PstateForOptimz = sensorData.pHistoryOracle(:,tNowIdx:sensorData.epochEndIdx);
    sensorData.PstateObjCalc = sensorData.pHistoryOracle(:,tNowIdx:sensorData.epochEndIdx);
end

% "Predict" fading coefficients
sensorData.commParam.h     = sensorData.commParam.hHistoryOracle(:,tNowIdx:sensorData.epochEndIdx);
sensorData.commParam.hactual= sensorData.commParam.hHistoryOracle(:,tNowIdx:sensorData.epochEndIdx);
% "Predict" energy usage
sensorData.energyHarvestForOptimz      = sensorData.eHistoryOracle(:,tNowIdx:sensorData.epochEndIdx);
sensorData.energyHarvestForOptimz(:,1) = sensorData.energyHarvestForOptimz(:,1) + sensorData.energyRem(:);
if any(any(sensorData.Emax < sensorData.energyHarvestForOptimz))
    sensorData.energyHarvestForOptimz = min(sensorData.Emax,sensorData.energyHarvestForOptimz);
    disp('WARNING: a least one starting energy value was too large and was reduced to Emax')
end

timeTemp = [time time(end)+(time(end)-time(end-1))];
sensorData.currentInterval = timeTemp(tNowIdx+1:1+sensorData.epochEndIdx) - timeTemp(tNowIdx:sensorData.epochEndIdx);