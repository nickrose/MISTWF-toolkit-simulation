function quantzVar = quantzNoiseFunc(bits, commParam,varargin)

% quantzNoiseFunc.m - given the number of bits and a set of communication
% parameters, returns the noise variance induced by quantization.
%
% USAGE:
%      quantzVar = quantzNoiseFunc(bits, commParam,varargin)
%   
% INPUTS:
%   bits                number of bits used for quantization
%   commParam           structure of communications parameters
%   varParamLengDiff    (optional) indicates whether the length of the bits
%                       variable is smaller than the dynamic range
%                       parameters by a multiple of d (dimension)
%    
% OUTPUTS:
%   quantzVar           the quantization variance
%
%
% Copyright (C) 2010  Nick Roseveare
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
% By: Nick Roseveare, 7-2010

if nargin < 3 || ~varargin{1}
    if length(bits) < length(commParam.W)
        numElem = length(commParam.W)/length(bits);
        bitsAll = repmat(bits',numElem,1); bitsAll = bitsAll(:);
        quantzVar = (commParam.W).^2 ./(3*(2.^bitsAll-1).^2);
    else
        quantzVar = (commParam.W).^2 ./(3*(2.^bits-1).^2);
    end
elseif varargin{1} > 1
    quantzVar = repmat(max(commParam.W),length(bits),1).^2 ./(3*(2.^bits-1).^2);
else
    numElem = length(commParam.W)/length(bits);
    indx = 1:numElem:length(commParam.W);
    quantzVar = (commParam.W(indx)).^2 ./(3*(2.^bits-1).^2);
end