function D = ObjFunc_eh_rwf(q)
% ObjFunc_eh.m - this objective function it meant to be used in conjunction
% with the energy harvesting optimization simulations.
% 
% USAGE:
%    D = ObjFunc_eh_rwf(q)
%
% Where 'q' is the input vector of optimization variables and 'Dinv' is the
% output distortion measure.
%
% The WF-WF objective
%   D = ( sum_k sum_n    q_nk /(P_n(k)*q_nk + W) )^(-1)
% The RWF-RWF objective
%   D =   sum_k sum_n    (P_n(k)*q_nk + W)/W^2
%
% Copyright (C) 2012  Nick Roseveare
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
% Last modified: January 2012

global sensorData

if length(q) == sensorData.numSensors*sensorData.numElem*sensorData.horizon
    % check if running single- or multi-dimension state estimate
    W = repmat(sensorData.commParam.W,sensorData.horizon,1);
else
    W = repmat(sensorData.commParam.Wp,sensorData.numSensors*sensorData.horizon,1);
end


D = sum((sensorData.PstateForOptimz(:).*q + W.^2)./q);
