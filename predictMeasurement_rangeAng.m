function [zhat, R, t] = predictMeasurement_rangeAng(sensorState,predictedState)
% predictMeasurement_rangeAng.m - Take the truth state and output a
% noisless prediction of the measurement given the input state.
% 
% USAGE:
% [zhat, R, t] = predictMeasurement_rangeAng(sensorState,predictedState)
%
% INPUTS:
%   sensorState  the sensor state structure which include measurement
%                function and covariance details (H and R)
%   predictedState      structure of truth trajectory information
%
% OUTPUTS:
%   zhat         the measurement vector
%   R            measurement covariance
%   t            measurement time
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
% By Nick Roseveare, 12-2011

H = sensorState.H;
measDim = size(H,1);

zhat = H*predictedState.xtrue;

sigR2  = sensorState.measNoiseVar.sigR2;
sigTh2 = sensorState.measNoiseVar.sigTh2;

phi = atan((zhat(2) - sensorState.senPos(2))/(zhat(1)-sensorState.senPos(1))) + (-0.5*sign(zhat(1)-sensorState.senPos(1))+0.5)*pi;
dist = norm(zhat-sensorState.senPos);
if sensorState.useRealisticRngAng
    R = @(r,theta)...
        [
        (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4)*(cos(theta))^2+(r^2*sigTh2+r*sigR2*sigTh2)*(sin(theta))^2,...
        (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4 - r^2*sigTh2-r*sigR2*sigTh2)*sin(theta)*cos(theta)...
        ;
        (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4 - r^2*sigTh2-r*sigR2*sigTh2)*sin(theta)*cos(theta)...
        (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4)*(sin(theta))^2+(r^2*sigTh2+r*sigR2*sigTh2)*(cos(theta))^2,...
        ];
    R = R(dist,phi);
else
    R = eye(measDim)*sensorState.measNoiseVar;
end

t = predictedState.time;
