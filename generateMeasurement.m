function [z, R, t] = generateMeasurement(H,truesys,sigma2v,tNowIdx,reuseNoise,runIdx)
% generateMeasurement.m - Take the truth state and output a noisy
% measurement. 
% 
% USAGE:
% [z, R, t] = generateMeasurement(H,truesys,sigma2v,tNowIdx,reuseNoise,mcRunIdx)
%
% INPUTS:
%   H            the observation function matrix
%   truesys      structure of truth trajectory information
%   sigma2v      variance of the measurement noise
%   tNowIdx      current time index to retrieve the current truth position
%                (note: not a time value, only an index)
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
% Copyright (C) 2009  Nick Roseveare
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

measDim = size(H,1);

% Database functions for reusing noise values
vSampleIdx = 1+(tNowIdx-1)*measDim:tNowIdx*measDim;
rvSample   = fixedSampleDatabase('measNoise',@randn,reuseNoise,vSampleIdx,runIdx);

v = rvSample'*sqrt(sigma2v);% sigma2v: measurement noise variance, 
z = H*truesys.xtrue(:,tNowIdx) + v; 
R = eye(measDim)*sigma2v;
t = truesys.tTrue(tNowIdx);