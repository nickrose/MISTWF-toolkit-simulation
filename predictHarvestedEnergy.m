function sensorData = predictHarvestedEnergy(scenario,sensorData,tNowIdx)
% predictHarvestedEnergy.m - provides a set of predicted
% amounts of energy available, and wireless network records the actual
% values in the time step following their realization. 
%
% USAGE:
%   sensorData = predictHarvestedEnergy(scenario,sensorData,currentTime)
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
% By: Nick Roseveare, October 2011

sensorData.energyHarvestForOptimz = [];
for nn = 1:sensorData.numSensors
    % Determine which MC state for energy quanta arrival this sensor is in
    MCstateVec = (scenario.energyHarvester.arrivalRateList-scenario.energyHarvester.arrivalRate(nn))==0 + 0;
    for kk = 1:sensorData.numEpochs
        MCstateVec = scenario.energyHarvester.arrivalRateMCtransition*MCstateVec;
        harvestedEnergyPred = MCstateVec'*scenario.energyHarvester.arrivalRateList * scenario.tsample*mean(scenario.energyHarvester.maxMinEnergy);
        sensorData.energyHarvestForOptimz(nn,kk) = harvestedEnergyPred;
    end
end

% Add the current battery levels to the first time step
firstEpochEnergy = sensorData.energyHarvestForOptimz(:,1);
sensorData.energyHarvestForOptimz(:,1) = sensorData.energyHarvestForOptimz(:,1) + sensorData.energyRem(:);
if any(any(sensorData.Emax < sensorData.energyHarvestForOptimz))
    sensorData.energyHarvestForOptimz = min(sensorData.Emax,sensorData.energyHarvestForOptimz);
    sensorData.energyRem = sensorData.energyRem - (firstEpochEnergy+sensorData.energyRem-sensorData.energyHarvestForOptimz(:,1));
    disp('WARNING: a least one predicted energy value was too large and was reduced to Emax')
end