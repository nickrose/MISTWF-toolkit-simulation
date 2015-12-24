function sensorData = predictFadingCoeff(scenario,sensorData,tIdx,time)
% predictFadingCoeff.m - assuming perfect CSI, this provides a
% set of predicted fading coefficients, and wireless network
% records the actual values in the time step following their realization.
%
% USAGE:
%   sensorData = predictFadingCoeff(scenario,sensorData,tIdx,time)
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

if isempty(sensorData.commParam.hHistory.time)
    % cannot make a prediction based on any information as none has been obtained
    sensorData.commParam.h = 1*ones(sensorData.numSensors,sensorData.numEpochs)+1e-6*randn(sensorData.numSensors,sensorData.numEpochs); % default value of avg power of Rayliegh fading factor ( E{|h_i|^2} )
    
else
    % Predict based on correlation what fading previously received
    timePast  = sensorData.commParam.hHistory.time;
    predSteps = time(tIdx:sensorData.epochEndIdx);
    
    Rxx = retrieveGaussianCorrFromRayl(timePast,timePast,sensorData.commParam.fd,'override');
    Ryx = retrieveGaussianCorrFromRayl(timePast,predSteps,sensorData.commParam.fd,'override');
    
    sensorData.commParam.h = [];
    for nn = 1:sensorData.numSensors
        sensorData.commParam.h(nn,:) = [(Ryx/Rxx)*sensorData.commParam.hHistory.prevH(nn,:)']';
        
    end
  
end
sensorData.commParam.hactual = sensorData.commParam.hHistoryOracle(:,tIdx:sensorData.epochEndIdx);


sensorData.commParam.hHistory.prevH = [sensorData.commParam.hHistoryOracle(:,tIdx) ... 
                            sensorData.commParam.hHistory.prevH(:,1:min(scenario.numSamplesPredBased-1,size(sensorData.commParam.hHistory.prevH,2)))...
                                        ];
sensorData.commParam.hHistory.time  = [time(tIdx) ...
                            sensorData.commParam.hHistory.time(1:min(scenario.numSamplesPredBased-1,size(sensorData.commParam.hHistory.time,2)))...
                                        ];
