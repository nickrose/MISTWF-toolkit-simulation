function [F,Q] = makeFQ(dim,tdiff,q,varargin)
% This function forms the forward state propagation matrix and the
% covariance process noise inflation matrix for a PV state with 'dim'
% dimensions and for the time difference 'tdiff'.
%
% Inputs:
%        dim  : the dimension of the PV state (1,2,3)
%        tdiff: the time difference from previous to current state
%        q    : the process noise variance
% Outputs:
%        F: state propagation matrix
%        Q: process noise inflation of covariance matrix
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
% By: Nick Roseveare, 12-2009

% State evolution matrix
if dim == 2
    if nargin == 3 || varargin{1} == 0% PV 
        F = [eye(dim)               tdiff*eye(dim);... % P components
            zeros(dim,dim)                 eye(dim)];   % V components
        
        Q = [1/3*tdiff^4*eye(dim), 1/2*tdiff^3*eye(dim);...
            1/2*tdiff^3*eye(dim),       tdiff^2*eye(dim)] * q;
        
    else % PVA
        F = [eye(dim)               tdiff*eye(dim)   tdiff^2*eye(dim);... % P components
            zeros(dim)              eye(dim)           tdiff*eye(dim);...
            zeros(dim)              zeros(dim)              eye(dim)];
        
        Q = [1/20*tdiff^5*eye(dim), 1/8*tdiff^4*eye(dim) 1/6*tdiff^3*eye(dim);...
            1/8*tdiff^4*eye(dim),  1/3*tdiff^3*eye(dim)     1/2*tdiff^2*eye(dim);...
            1/6*tdiff^3*eye(dim)   1/2*tdiff^2*eye(dim)      tdiff*eye(dim)] * q;
    end
end