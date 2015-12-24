function sensorData = predictAllEstimate(sensorData,sensorState,scenario,tNowIdx,time)
% predictAll.m
%
% USAGE:
%     sensorData = predictAllEstimate(sensorData,sensorState,scenario,tNowIdx,time)
%
% sensorData        the parameter and variable structure holding most of
%                   the operational data for the system and esp. for the
%                   optimization 
% sensorState       holds only the current state model and history of state
%                   updates and measurements as well as their covariances
% tNowIdx           a unique index into the scenario simulation time vector
% time              a vector of values nsamples long with tsample step
%                   spacing between samples
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

% Predict state and sensor noise variances

sensorData.PstateForOptimz = [];
if scenario.usesObjectTraj
    sensorData = predictStateAndSensorUncertainty(sensorData,sensorState,scenario,tNowIdx,time);
end
% Predict fading coefficients
sensorData = predictFadingCoeff(scenario,sensorData,tNowIdx,time);

% Predict energy usage
sensorData = predictHarvestedEnergy(scenario,sensorData,tNowIdx);

timeTemp = [time time(end)+(time(end)-time(end-1))];
sensorData.currentInterval = timeTemp(tNowIdx+1:1+min(tNowIdx+sensorData.horizon-1,length(time))) - timeTemp(tNowIdx:sensorData.epochEndIdx);