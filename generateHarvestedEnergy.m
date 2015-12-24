function sensorData = generateHarvestedEnergy(scenario,sensorData,time)
% generateHarvestedEnergy.m - assuming perfect battery
% information is fed back, this provides a set of actual 
% amounts of energy available, and wireless network records the actual
% values in the time step following their realization. 
%
% USAGE:
%   sensorData = generateHarvestedEnergy(scenario,sensorData,time)
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
global recycleMCdraws


arrivalRate = scenario.energyHarvester.arrivalRate;
prevTime = -scenario.tsample;
timeLength = length(time);
for kk = 1:timeLength
    if sum(diag(scenario.energyHarvester.arrivalRateMCtransition)) ~= length(scenario.energyHarvester.arrivalRateList)
    % This condition avoids transitioning arrival rates if the transition matrix is the identity matrix
        
        for nn = 1:scenario.numSensors
            % Determine which MC state for energy quanta arrival this sensor is in
            MCstateVec = abs(scenario.energyHarvester.arrivalRateList-arrivalRate(nn))<1e-10 + 0;
            
            if kk>1
                % Determine the Markov chain transitions
                threshold = fixedSampleDatabase('markovChainTrans',@rand,...
                    recycleMCdraws.reuseNoise,kk+(nn-1)*timeLength,...
                    recycleMCdraws.runIdx);
                transProb = MCstateVec'*scenario.energyHarvester.arrivalRateMCtransition;
                MCstateVec = zeros(length(MCstateVec),1);
                for pp = 1:length(MCstateVec)
                    if threshold <= sum(transProb(1:pp))
                        MCstateVec(pp) = 1;
                        break
                    end
                end
                arrivalRate(nn) = scenario.energyHarvester.arrivalRateList(MCstateVec==1);
            end
        end
    end
    if kk>1
        prevTime = time(kk-1);
    end
    % Harvest the energy for this time step
    randEnergyLevels = fixedSampleDatabase('energyLevels',@rand,...
                    recycleMCdraws.reuseNoise,1+(kk-1)*scenario.numSensors:kk*scenario.numSensors,...
                    recycleMCdraws.runIdx)';
                
    poissonVec = poissrnd(arrivalRate*(time(kk)-prevTime));
    arrivalRatesRVs = fixedSampleDatabase('energyArrivals','data',...
                    recycleMCdraws.reuseNoise,1+(kk-1)*scenario.numSensors:kk*scenario.numSensors,...
                    recycleMCdraws.runIdx,...
                    poissonVec);

    
    sensorData.eHistoryOracle(:,kk) = arrivalRatesRVs'.*...
        ([ones(scenario.numSensors,1), randEnergyLevels]*...
            [scenario.energyHarvester.maxMinEnergy(1);...
             scenario.energyHarvester.maxMinEnergy(2)-scenario.energyHarvester.maxMinEnergy(1)]);
end

