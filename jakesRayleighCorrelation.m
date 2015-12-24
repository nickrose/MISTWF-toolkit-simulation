function [corr,varargout] = jakesRayleighCorrelation(fd,time,varargin)
% jakesRayleighCorrelation.m - find correlation between samples instances
% accorading to Jakes' model.
%
% NOTE: 
% If the product of the max doppler frequency and the time step are in
% accordance with the following ranges, then the expected behaviour of the
% correlation is 
%             2*pi*fd*timeStep < .01      Highly correlated (~1)
%         2 < 2*pi*fd*timeStep < 10       Less correlated (periodic, pos & neg)
%        10 < 2*pi*fd*timeStep            Nearly uncorrelated (~0)
%
% Alternatively it is useful to put the constants at the limits
%             fd*timeStep < 0.0015      Highly correlated (~1)
%     0.318 < fd*timeStep < 1.592       Less correlated (periodic, pos & neg)
%     1.592 < fd*timeStep               Nearly uncorrelated (~0)

%             
% USAGE:
%   [corr,tau] = jakesRayleighCorrelation(fd,time)
%   [corr,tau] = jakesRayleighCorrelation(fd,time,nsamples)
%    corr      = jakesRayleighCorrelation(fd,timeStepIn,timeStepOut)
%    corr      = jakesRayleighCorrelation(fd,timeStepIn,timeStepOut,cmd)
%
% INPUTS:
%   fd              Max Doppler frequency (Hertz) (1/fd is coh time)
%   time            The time difference scalar or vector of time
%                   differences (assumes a single realized values and
%                   computes standard correltion function) (seconds)
%   timeStepIn      Can be either a scalar representing a uniform time
%   timeStepOut     difference, or a vector of time instances;
%                   In and Out represent the number of realized samples and
%                   the samples being predicted, respectively
%   nsamples        (optional) if timeStep is a scalar then return a set of
%                   equispaced correlation values
%   cmd             (optional) can be 
%                   'override'  for the case when timeStepIn and timeStep
%                   out are scalars, but nsamples is NOT a desired input,
%                   or
%                   'useApprox' approximation where appropriate, (T/F)
%
% OUTPUTS:
%   corr        the correlation values (scaled between 0 and 1)
%   tau         (optional) the generated vector equispace time steps
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
% By: Nick Roseveare, October, 2011

if nargin > 2 
    if ~isscalar(varargin{1}) || (~isempty(varargin{2}) && strcmp(varargin{2},'override'))
        timeStepIn = time;
        timeStepOut = varargin{1};
        tau = abs(repmat(timeStepOut',1,length(timeStepIn)) - repmat(timeStepIn,length(timeStepOut),1)); % abs value represents the fact that only the amount of time different matters
    else
        nsamples = varargin{1};
        if ~isscalar(time)
            error('time variable should not be a vector when specifying number of samples')
        end
        tau = (1:(nsamples-1))*time;
    end
else
    tau = abs(time);
end

% Approximations for possible future speed ups
if nargin == 4 && (~isempty(varargin{2}) && strcmp(varargin{2},'useApprox')) % <- useApprox
    tauColumize = tau(:);
    corr = zeros(length(tauColumize),1);
    for kk = 1:length(tauColumize)
        if  2*pi*fd*tauColumize(kk) < .01
            corr(kk) = 1;
        elseif 2*pi*fd*tauColumize(kk) > 2
            corr(kk) = sqrt(2/(pi*2*pi*fd*tauColumize(kk)))*cos(2*pi*fd*tauColumize(kk) - pi/4);
        else
            corr(kk) = besselj(0,2*pi*fd*tauColumize(kk));
        end
    end
    corr = reshape(corr,size(tau));
else
    % Just use matlab built-in bessel function
    corr = besselj(0,2*pi*fd*tau);
end

if length(time) == 1 && ~(nargin > 3 && (~isempty(varargin{2}) && strcmp(varargin{2},'override')))
    if nargout == 2
        varargout{1} = tau;
    else
        varargout{1} = [];
    end
end
