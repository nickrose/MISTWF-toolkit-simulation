function [R, corrSign] = retrieveGaussianCorrFromRayl(timeStepIn,timeStepOut,coherencebw,varargin)
% retrieveGaussianCorrFromRayl.m
%
% This function determines the cross corelation between time instances and
% then calculates the corresponding correlation between the complex gaussian 
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
% By: Nick Roseveare, Nov 2011
% Modified: Oct 2012 - error (NaN from interpolation) resulted since values
% of gaussian correlation previously only went down to 0.0018, changed the
% previous value of zero to something 'visually' interpolated' between the
% next value up and zero. See previous list of values here:
%
% x = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1];
% y = [0.0018 0.0047 0.0056 0.0243 0.0337 0.0559 0.0737 0.0965 0.1494 0.1836 0.2227 0.2752 0.3327 0.4133 0.4562 0.5410 0.6073 0.6974 0.7913 0.9005 1]; 
%
% Modeled after 'correl.m' by B. Natarajan

x = [0 0.028 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1];
y = [0 0.0018 0.0047 0.0056 0.0243 0.0337 0.0559 0.0737 0.0965 0.1494 0.1836 0.2227 0.2752 0.3327 0.4133 0.4562 0.5410 0.6073 0.6974 0.7913 0.9005 1]; 

% /\ addded a (0,0) to vectors above to prevent NaN resulting after
% interpolation


%rho = 1/(1+((i-1)*deltaf/coherencebw)^2) - pi/(4*N);
if nargin > 3
    rho = jakesRayleighCorrelation(coherencebw,timeStepIn,timeStepOut,varargin{1});
else
    rho = jakesRayleighCorrelation(coherencebw,timeStepIn,timeStepOut);
end
[m,n] = size(rho);

rho = rho(:);

rhonew = zeros(1,length(rho));
corrSign = zeros(size(rhonew));
for i = 1:length(rho),
    % Attempting to preserve the sign, not sure if this is correct
    rhonew(i) = interp1(y,x,abs(rho(i))); 
    corrSign(i) = sign(rho(i));
end

R = reshape(rhonew,m,n); % return the normalized correlation matrix
corrSign = reshape(corrSign,m,n);
