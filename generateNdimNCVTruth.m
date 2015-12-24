function [truesys, H, G] = generateNdimNCVTruth(dim,initState,q,tsample,nsamplesIN,obsRegion)
% This function takes in the sample rate and the number of samples and
% generates the system state space truth data for a 1-dim constant velocity 
% motion model.
% 
% Inputs: 
%         dim        : number of dimensions for PV state (1,2,3)
%         initState  : the initial state 
%         q          : the spectral density of the process noise (variance)
%         tsample    : time indices OR sample rate
%         nsamplesIN : if 'tsample' is a rate, then must specify the number
%                      of samples to generate truth over 
%         obsRegion  : (optional) the region in which the generated truth
%                      path should stay, the fifth element is the max
%                      magnitude of the velocity
%
% Outputs:
%         truesys  : a structure variable with sub fields specifying the
%                    true dynamics of the generated sys
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

obsRegionUsed = nargin == 6;

numElemPV = dim*2; % number of elements in the PV for this number of dimensions

G =  eye(numElemPV); % 'B' matrix in standard SS, identity implies independent noise in the position and velocity components
H =  [eye(dim) zeros(dim)];  % observation matrix (observes position only)
x0 = initState; % true initial state

x = x0;
P = blkdiag(eye(dim),eye(dim)/500)*q;% using this is CRLB, change there also



if obsRegionUsed % generate truth within a constrained region
    nsamples = nsamplesIN;
    
    tTrue = 0:tsample:(nsamples-1)*tsample;
    
    % Only need to form the F and Q matrices once
    [F,Q] = makeFQ(dim,tsample,q);
    
    restart = 1;
    outerCounter = 0;
    while restart
        outerCounter = outerCounter + 1;
        [xtrue,Ptrue,restart] = createStateEvolBndd(x,P,F,Q,G,obsRegion,numElemPV,nsamples);
    end
    
    % Form the true system struct
    truesys.Ptrue = Ptrue;
    truesys.xtrue = xtrue;
    truesys.tTrue = tTrue;
    
elseif nargin <= 5  % unbounded generation of truth data
    if nargin == 5
        nsamples = nsamplesIN;
        if length(tsample) > 1
            warning('Do not need to provide a sample rate in addition to a list of sample indices\n')
        end
    else
        nsamples = length(tsample);
    end
    % Test to see if the samples vector is uniformly spaced
    if all((tsample(2:end)-tsample(1:end-1) - tsample(3:end)-tsample(2:end-1)) < 1e-8);
        % Call self with uniformly spaced sampling to speed up
        [truesys, H, G] = generateNdimNCVTruth(dim,initState,q,tsample,nsamples);
    else% user determined sample spacing
        
        % Proceed normally with non-uniformly spaced samples
        [xtrue,Ptrue] = createStateEvol(x,P,q,G,numElemPV,nsamples,tsample);
        
        % Form the true system struct
        truesys.Ptrue = Ptrue;
        truesys.xtrue = xtrue;
        truesys.tTrue = tsample;
    end
else
    error('Incorrect number of inputs')
end
% ===============================================================

function [xtrue,Ptrue,restart] = createStateEvolBndd(x,P,F,Q,G,obsRegion,numElemPV,nsamples)
% Create the evolving system states.
global recycleMCdraws

xtrue = zeros(numElemPV,nsamples);
Ptrue = zeros(numElemPV,numElemPV,nsamples);

% idx = 1:numElemPV;
% procNoise = randn(1,numElemPV);
% rvDraw = fixedSampleDatabase('trajProcNoiseBdd','data',recycleMCdraws.reuseNoise,idx,recycleMCdraws.runIdx,procNoise);
% useDatabase = (norm(procNoise-rvDraw)>1e-12);

for n = 1:nsamples
    % Record the current state
    xtrue(:,n) = x;
    Ptrue(:,:,n) = P;
    
    newPointOkay = 0;
    xPrev = x;
    % Propagate the state
    Q1half = sqrtm(Q);
    
%     idx = (1+(n-1)*numElemPV):n*numElemPV;
    tryCount = 0;
    while ~newPointOkay && tryCount < 50
        procNoise = randn(1,numElemPV);
%         if useDatabase
%             rvDraw = fixedSampleDatabase('trajProcNoiseBdd','data',recycleMCdraws.reuseNoise,idx,recycleMCdraws.runIdx,procNoise);
%         else
            rvDraw = procNoise;
%         end
        
        w = Q1half*rvDraw';  % generate process noise
        x = F*xPrev + G*w; % propagate the state
        
        
        % NOTE: works only for dim = 2 right now, would need to modify
        % for other dimensions
        newPointOkay = x(1) >= obsRegion(1) && ...
            x(1) < obsRegion(2) && ...
            x(2) >= obsRegion(3) && ...
            x(2) < obsRegion(4) &&...
            norm(x(3:4)) < obsRegion(5); % ensure velocity is constrained
        tryCount = tryCount + 1;
    end
%     if ~useDatabase
%         % record data into database - didn't want to waste time waiting for
%         % an appropriate rv sample in the while loop
%         fixedSampleDatabase('trajProcNoiseBdd','data',recycleMCdraws.reuseNoise,idx,recycleMCdraws.runIdx,procNoise);
%     end
    P = Q;%F*P*F' +  Q  % don't propagate the covariance, this is the truth, so it should have the same covariance at every point
    
    restart = ~newPointOkay;
    if ~newPointOkay % don't want to run an infinite loop, just restart if the process won't cooperate
        break;
    end
end


% ===============================================================

function [xtrue,Ptrue] = createStateEvol(x,P,q,G,numElemPV,nsamples,tsample)
% Create the evolving system states.
global recycleMCdraws

xtrue = zeros(numElemPV,nsamples);
Ptrue = zeros(numElemPV,numElemPV,nsamples);

% Record the current state
xtrue(:,1) = x;
Ptrue(:,:,1) = P;

% Create the evolving system states.
for n = 2:nsamples
    
    % Need to form the F and Q matrices each iteration because time spacing
    % not uniform
    tdiff = tsample(n) - tsample(n-1);
    [F,Q] = makeFQ(dim,tdiff,q);
    
    xPrev = x;
    
    idx = (1+(n-1)*numElemPV):n*numElemPV;
    rvDraw = fixedSampleDatabase('trajProcNoiseUnbdd',@randn,recycleMCdraws.reuseNoise,idx,recycleMCdraws.runIdx);
    
    
    % Propagate the state
    w = sqrtm(Q)*rvDraw';  % generate process noise
    x = F*xPrev + G*w; % propagate the state
    
    % Record the current state
    xtrue(:,n) = x;
    Ptrue(:,:,n) = P;
    
end

