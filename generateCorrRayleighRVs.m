function alphas = generateCorrRayleighRVs(R,corrSign,idxStr)
% generateCorrRayleighRVs.m
%
% This function takes the required NxN covariance matrix (R) of the complex
% gaussian RVs and generates N correlated Rayleigh random variables. 
%
% Variance of gaussian = 1;
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
% By: Nick Roseveare, Dec 2011
% Modeled after 'raycor2.m' by B. Natarajan
global recycleMCdraws

if nargin < 3
    idxStr.curIdx = 1;
    idxStr.setLen = 1;
end

N = size(R,1);
xSampleIdx = (1+(idxStr.curIdx-1)*idxStr.setLen):(idxStr.curIdx*idxStr.setLen);
ySampleIdx = (1+(idxStr.curIdx-1)*idxStr.setLen):(idxStr.curIdx*idxStr.setLen);
rvSampleX = fixedSampleDatabase('corrRayNoiseI',@randn,recycleMCdraws.reuseNoise,xSampleIdx,recycleMCdraws.runIdx);
rvSampleY = fixedSampleDatabase('corrRayNoiseQ',@randn,recycleMCdraws.reuseNoise,ySampleIdx,recycleMCdraws.runIdx);

% First Generate K uncorrelated gaussian RVs
x = rvSampleX'*sqrt(0.5);
y = rvSampleY'*sqrt(0.5);
z = x + 1i*y;

% Perform Cholesky Decomposition
    [U,Z] = eig(R);
    Z = diag(sqrt(diag(Z)));
    D = Z*U';
    % T = inv(D);
% L = chol(R,'upper'); % cholesky factorize because matrix Pos Semi def, cheaper "inverse" (i.e. back solve)


% Transformation
cmplxEnv = D\z;% T*z
% cmplxEnv = L\z; % essentially a 'unwhitening' procedure y = R^(-1/2)*x

% Find attenuation envelope
alphas = abs(cmplxEnv);

