function U = ObjFunc_eh_wtlog(p)
% ObjFunc_eh_wtlog.m - this objective function it meant to be used in
% conjunction with the energy harvesting optimization simulations.
% 
% USAGE:
%    U = ObjFunc_eh_wtlog(p)
%
% Where 'p' is the input vector of optimization variables and 'U' is the
% sum throughput rate utility function.
%
% The objective
%   U = sum_k sum_n    w_nk * log(1 + h_nk * p_nk)
%
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
% By: Nick Roseveare, March 2012
% Borrowed from:
% ObjFunc_eh.m, by: Nick Roseveare, December 2011

global sensorData

U = -sum(sensorData.weights(:).*log(1+sensorData.commParam.h(:).*p));