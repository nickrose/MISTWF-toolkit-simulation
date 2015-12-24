function [z, R, t] = generateMeasurement_rangeAng(sensorState,truesys,tNowIdx,useRealisticRngAng)
% generateMeasurement_rangeAng.m - Take the truth state and output a noisy
% measurement according with a correlated range angle covariance. 
% 
% USAGE:
% [z, R, t] = generateMeasurement_rangeAng(H,truesys,tNowIdx,useRealisticRngAng)
%
% INPUTS:
%   sensorState  the sensor state structure which include measurement
%                function and covariance details (H and R)
%   truesys      structure of truth trajectory information
%   tNowIdx      current time index to retrieve the current truth position
%                (note: not a time value, only an index)
%   useRealisticRngAng   speaks for itself
%<not used, uses global variable instead>
%   reuseNoise   used by the simulation to reuse random numbers across
%                simulation instances
%   runIdx       Monte Carlo/scenario run index, used as a part of the previous
%
% OUTPUTS:
%   z            the measurement vector
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
% By Nick Roseveare, 12-2009
% Updated: 2-2011

global recycleMCdraws

H = sensorState.H;
measDim = size(H,1);

% Database functions for reusing noise values
vSampleIdx = 1+(tNowIdx-1)*measDim:tNowIdx*measDim;
if ~isempty(recycleMCdraws)
    rvSample = fixedSampleDatabase('measNoise',@randn,recycleMCdraws.reuseNoise,vSampleIdx,recycleMCdraws.runIdx);
    rvSample = rvSample';
else
    rvSample = randn(measDim,1);
end
if numel(truesys.xtrue) < size(H,2)
    z = truesys.xtrue(:,tNowIdx);
else
    z = H*truesys.xtrue(:,tNowIdx);
end

sigR2  = sensorState.measNoiseVar.sigR2;
sigTh2 = sensorState.measNoiseVar.sigTh2;

phi = atan((z(2) - sensorState.senPos(2))/(z(1)-sensorState.senPos(1))) + (-0.5*sign(z(1)-sensorState.senPos(1))+0.5)*pi;
dist = norm(z-sensorState.senPos);
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
    v = sqrtm(R)*rvSample;% measurement noise covariance, 
else
    R = eye(measDim)*sensorState.measNoiseVar;
    v = sqrtm(R)*rvSample;
end

z = z + v; 
t = truesys.tTrue(tNowIdx);