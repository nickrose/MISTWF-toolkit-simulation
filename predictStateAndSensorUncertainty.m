function sensorData = predictStateAndSensorUncertainty(sensorData,sensorState,scenario,tNowIdx,time)
% predictStateAndSensorUncertainty.m - predict the uncertainty associated
% with the measurements and state update for future time instances. 
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

predictedSensorState = sensorState;
sensorData.PstateObjCalc = sensorData.pHistoryOracle(:,tNowIdx:sensorData.epochEndIdx);

sensorData.PstateForOptimz = [];
for jj = tNowIdx:sensorData.epochEndIdx 
    for nn = 1:sensorData.numSensors
        if jj > tNowIdx
            [F,Q] = makeFQ(scenario.numDimen,time(jj)-predictedSensorState(nn).t,predictedSensorState(nn).q);
            if time(jj)-predictedSensorState(nn).t > 1e-8
                % Update state to current time
                predictMeas.xtrue = F*predictedSensorState(nn).x;       % x[n+1|n]
                predictMeas.time = time(jj);
            end
            [z, R, t] = predictMeasurement_rangeAng(predictedSensorState(nn),predictMeas);
            predictedSensorState(nn) = filterMeasurement(predictedSensorState(nn),z,R,t);
        end
        if sensorData.pMetricSize == 1
            if jj == tNowIdx
                sensorData.PstateForOptimz(nn,jj-tNowIdx+1) = scenario.optFormSetup.PmetricFunc(predictedSensorState(nn).P);
            else
                % there is a j "+ 1" because the "initial" sensor state is
                % before the "0 state" - updated after the 0th measurement
                sensorData.PstateForOptimz(nn,jj-tNowIdx+1) = scenario.optFormSetup.PmetricFunc(predictedSensorState(nn).PHist(:,:,jj+1));
            end
        else
            sensorData.PstateForOptimz(1+numElem*(nn-1):(numElem+nn),jj-tNowIdx+1) = scenario.optFormSetup.PmetricFunc(predictedSensorState(nn).P); % state covariance history
        end
    end
end