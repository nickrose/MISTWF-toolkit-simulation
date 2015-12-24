function sensorState = makeSensorStructEH(scenario,time,F,Q)
% Creates a vector struct of sensor-level Kalman filter information.
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
global recycleMCdraws

numSensors = scenario.numSensors;
dim     = scenario.numDimen;
H       = scenario.trueSys.H;
q       = scenario.q;
p0      = scenario.p0;

numElements = dim*2;
sensorState = struct('senPos',zeros(dim,1),'distToFC',[],'P',zeros(numElements,numElements),'x',zeros(numElements,1),'t',[],'H',[],'measNoiseCovFunc',[]);
sensorState = repmat(sensorState,numSensors,1);
Hsize = size(H,1);


if length(time) > 1
    tnow = time(1);
else
    tnow = 0;
end
if isscalar(p0)
    p0 = p0*[eye(dim) zeros(dim);zeros(dim) eye(dim)*2];
end

newlyGenerated = 0;

% Make or retrieve the sensor positions
if scenario.useLastWSNpos && ~exist('lastSenPos','var') && exist([scenario.simDataSaveDirectory, 'lastTrueWSNpos.mat'],'file')
    load([scenario.simDataSaveDirectory, '/lastTrueWSNpos'])
end

if ~scenario.newWSNposEachMCrun && scenario.useLastWSNpos ...
        && exist('lastSenPos','var')...
        && ...
          (  ((size(lastSenPos(1).sensorPos,2) == scenario.numSensors) ...
              && (isempty(recycleMCdraws) || ~recycleMCdraws.reuseNoise))...
            || ...
             ((~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise)...
              && (length(lastSenPos) >= recycleMCdraws.runIdx))...
              && (size(lastSenPos(recycleMCdraws.runIdx).sensorPos,2) == scenario.numSensors)...
          )
    % check for the reusage of WSN positions
    if ~isempty(recycleMCdraws)&& recycleMCdraws.reuseNoise &&...
            (length(lastSenPos) < recycleMCdraws.runIdx)
        
        sensorPos = getSensorPos(scenario,numSensors);
        lastSenPos(recycleMCdraws.runIdx).sensorPos = sensorPos;
        
    elseif ~isempty(recycleMCdraws) &&...
            recycleMCdraws.reuseNoise
        sensorPos = lastSenPos(recycleMCdraws.runIdx).sensorPos;
    else
        sensorPos = lastSenPos(1).sensorPos;
    end
    
    save([scenario.simDataSaveDirectory, '/lastTrueWSNpos'],'lastSenPos')
    
elseif scenario.newWSNposEachMCrun && (isempty(recycleMCdraws) || (~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise))
    sensorPos = getSensorPos(scenario,numSensors);
    lastSenPosTemp.sensorPos = sensorPos;
    if ~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise
        lastSenPos(recycleMCdraws.runIdx) = lastSenPosTemp;
    else
        lastSenPos = lastSenPosTemp;
    end
    save([scenario.simDataSaveDirectory, '/lastTrueWSNpos'],'lastSenPos')
else
    sensorPos = getSensorPos(scenario,numSensors);
    newlyGenerated = 1;
    lastSenPosTemp.sensorPos = sensorPos;

    if ~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise
        lastSenPos(recycleMCdraws.runIdx) = lastSenPosTemp;
    else
        lastSenPos = lastSenPosTemp;
    end
    save([scenario.simDataSaveDirectory, '/lastTrueWSNpos'],'lastSenPos')
    
end

if ~scenario.useLastWSNpos || newlyGenerated
    figure(13);
    hold on;
    plot(sensorPos(1,:),sensorPos(2,:),'rx','MarkerSize',10)
    if scenario.useFirstNew
        title('Sensor positions');
    else
        title('Possible sensor positions okay? (zero (0) to cancel)');
    end
    axis(scenario.observRegion(1:4))
    axis equal
    hold off;
    if ~scenario.useFirstNew
    cont = input('enter 0 to cancel:  ');
        if cont ==0
            keyboard;
        end
    end
    pause(.5)
    save([scenario.simDataSaveDirectory, 'lastTrueWSNpos'],'lastSenPos')
elseif scenario.plotExtraSetupInfo
    figure(13);
    hold on;
    plot(sensorPos(1,:),sensorPos(2,:),'rx','MarkerSize',10);hold on;
    offset = scenario.observRegion([2 4])-scenario.observRegion([1 3]);
    for ii = 1:size(sensorPos,2)
        text(sensorPos(1,ii)+offset(1)/40,sensorPos(2,ii)+offset(2)/40,num2str(ii));hold on;
    end
    title('Sensor positions and Object trajectory');
    axis(scenario.observRegion(1:4))
    axis equal
    hold off;
    pause(.5)
    
end
n = scenario.noiseVar; %(normally 1, smaller for "nice" plots) base measurement noise variance

for kk = 1:numSensors
    sensorState(kk).x = zeros(numElements,1);
    sensorState(kk).P = p0;
    sensorState(kk).xHist = zeros(numElements,1);
    sensorState(kk).zHist = zeros(Hsize,1);
    sensorState(kk).PHist = [eye(dim) zeros(dim);zeros(dim) eye(dim)/10];
    sensorState(kk).RHist = eye(Hsize);
    sensorState(kk).tHist = tnow;
    sensorState(kk).H = H;
    sensorState(kk).q = q;
    sensorState(kk).t = tnow;
    sensorState(kk).senPos = sensorPos(:,kk);
    sensorState(kk).distToFC = norm(scenario.FCnodePosition-sensorPos(:,kk));
    sensorState(kk).useRealisticRngAng = scenario.useRealisticRngAng;
    
    if scenario.useRealisticRngAng
        sigR2  = scenario.rangeNoiseVar * scenario.noiseVar(kk);
        sensorState(kk).measNoiseVar.sigR2 = sigR2;
        sigTh2 = scenario.angNoiseVar * scenario.noiseVar(kk);
        sensorState(kk).measNoiseVar.sigTh2 = sigTh2;
        
        sensorState(kk).measNoiseCovFunc = @(r,theta)...
            [
            (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4)*(cos(theta))^2+(r^2*sigTh2+r*sigR2*sigTh2)*(sin(theta))^2,...
            (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4 - r^2*sigTh2-r*sigR2*sigTh2)*sin(theta)*cos(theta)...
            ;
            (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4 - r^2*sigTh2-r*sigR2*sigTh2)*sin(theta)*cos(theta)...
            (r*sigR2+r^2*(sin(sqrt(sigTh2)))^4)*(sin(theta))^2+(r^2*sigTh2+r*sigR2*sigTh2)*(cos(theta))^2,...
            ];
    else
        sensorState(kk).measNoiseVar = n;
        sensorState(kk).measNoiseCovFunc = [n(kk) 0;0 n(kk)];
    end
end
if nargin > 2
    for kk = 1:numSensors
        sensorState(kk).F = F;
        sensorState(kk).Q = Q;
    end
end

% ===========================================================
function sensorPos = getSensorPos(scenario,numSensors)

global recycleMCdraws

if ~isempty(recycleMCdraws)
    xVecPos = 1:numSensors;
    yVecPos = (1:numSensors) + numSensors;
    rvSampleX = fixedSampleDatabase('WSNpositions',@rand,recycleMCdraws.reuseNoise,xVecPos,recycleMCdraws.runIdx);
    rvSampleY = fixedSampleDatabase('WSNpositions',@rand,recycleMCdraws.reuseNoise,yVecPos,recycleMCdraws.runIdx);
else
    rvSampleX = rand(1,numSensors);
    rvSampleY = rand(1,numSensors);
end
sensorPos = zeros(2,numSensors);

if scenario.randDistWSNpos
    sensorPos(1,:) = rvSampleX*scenario.observRegion(2)+scenario.observRegion(1);
    sensorPos(2,:) = rvSampleY*scenario.observRegion(4)+scenario.observRegion(3);
else
    rowCol = factor(floor((numSensors+1)/2)*2);
    while length(rowCol) > 2
        factors = rowCol;
        rowCol = [];
        while ~isempty(factors)
            [newMin,idxm] = min(factors);
            [newMax,idxM] = min(factors);
            factors([idxm idxM])= [];
            rowCol = [rowCol newMin*newMax];
        end
    end
    for sensorsPlaced = 1:numSensors
        sensorPos(2,sensorsPlaced) = (1/(rowCol(1)+1)+ mod(sensorsPlaced-1,rowCol(1))/(rowCol(1)+1)+(2*rand(1,1)-1)/(rowCol(1)+1))*scenario.observRegion(4)+scenario.observRegion(3);
        sensorPos(1,sensorsPlaced) = (1/(rowCol(2)+1)*(1+floor(sensorsPlaced/(rowCol(1)+1))+(2*rand(1,1)-1)/(rowCol(2)+1)))*scenario.observRegion(2)+scenario.observRegion(1);
    end
end