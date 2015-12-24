function ploss = calcPathLoss(dist,PLexponent)
% calcPathLoss.m - simple pathloss calculation
%
% INPUTS:
%   dist        vector of distances from a fixed point (could be a matrix
%               with zeros along the diagonals
%   PLexponent  the path loss exponent (usually between 2 and 4 for some
%               what unreflective environment up through flat earth
%               reflection models, 4-6 for highly cluttered environments)
% OUTPUT:
%   ploss       the path loss (usually a small positive quantity, which
%               scales the transmission power, thus attenuating it)
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

ploss = dist.^(-PLexponent);