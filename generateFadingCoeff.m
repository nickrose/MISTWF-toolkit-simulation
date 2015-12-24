function sensorData = generateFadingCoeff(scenario,sensorData,currentTime)
% generateFadingCoeff.m - assuming perfect CSI, this provides a
% set of actual fading coefficients, and wireless network
% records the actual values in the time step following their realization.
%
% USAGE:
%   sensorData = generateFadingCoeff(scenario,sensorData,currentTime)
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
% By: Nick Roseveare, October 2011, not being used in current simul form.

if isempty(sensorData.hHistory.time)
    % cannot make a prediction based on any information as none has been obtained
    sensorData.commParam.hPredict = 1*ones(numSensors,sensorData.horizon); % default value of avg power of Rayliegh fading factor ( E{|h_i|^2} )
    
    % ===============================================================
    % Generate true fading
    sensorData.actualH.i = randn(sensorData.numSenors,1);
    sensorData.actualH.q = randn(sensorData.numSenors,1);
    sensorData.actualH.mag = sqrt(sensorData.actualH.i^2+sensorData.actualH.q^2);
else
  
    
    % ===============================================================
    % Generate true fading
    nI = randn(sensorData.numSenors,1);
    nQ = randn(sensorData.numSenors,1);
    corr = jakesRayleighCorrelation(sensorData.commParam.fd,[currentTime time],currentTime);
    for nn = 1:scenario.numSensors
        nI(nn) = corr*[nI(nn); prevI'];
        nQ(nn) = corr*[nQ(nn); prevQ'];
    end
    sensorData.actualH.i = nI;
    sensorData.actualH.q = nQ;
    sensorData.actualH.mag = sqrt(sensorData.actualH.i^2+sensorData.actualH.q^2);
end