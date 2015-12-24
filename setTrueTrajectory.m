function trueSys = setTrueTrajectory(scenario)
% setTrueTrajectory.m - make or retrieve true object trajectory
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
% By: Nick Roseveare, Oct 2011

global recycleMCdraws

if scenario.useLastTruth && ~exist('lastTrue','var') && exist([scenario.simDataSaveDirectory, '/lastTrueFile.mat'],'file')
    load([scenario.simDataSaveDirectory '/lastTrueFile'])
end

if ~scenario.newTrajEachMCrun && scenario.useLastTruth && exist('lastTrue','var') ...
        &&( ( ( isempty(recycleMCdraws) || ~recycleMCdraws.reuseNoise) ...
            && (lastTrue(1).qtrue    == scenario.qtrue) ...
            && (lastTrue(1).nsamples == scenario.nsamples) ...
            && (lastTrue(1).tsample  == scenario.tsample)...
            && (lastTrue(1).numDimen == scenario.numDimen)...
            && (lastTrue(1).numSensors == scenario.numSensors)...
            && (norm(lastTrue(1).initState - scenario.initState) < 1e-12)...
            )...
            ||...
            (  (~isempty(recycleMCdraws))...
            && (length(lastTrue) >= recycleMCdraws.runIdx)...
			&& ~isempty(lastTrue(recycleMCdraws.runIdx).nsamples)...
            && (lastTrue(recycleMCdraws.runIdx).qtrue    == scenario.qtrue) ...
            && (lastTrue(recycleMCdraws.runIdx).nsamples == scenario.nsamples) ...
            && (lastTrue(recycleMCdraws.runIdx).tsample  == scenario.tsample)...
            && (lastTrue(recycleMCdraws.runIdx).numDimen == scenario.numDimen)...
            && (lastTrue(recycleMCdraws.runIdx).numSensors == scenario.numSensors)...
            && (norm(lastTrue(1).initState - scenario.initState) < 1e-12)...
            )...
           )
    if ~isempty(recycleMCdraws)&& recycleMCdraws.reuseNoise &&...
            (length(lastTrue) < recycleMCdraws.runIdx)
        [truesys, H, G]     = generateNdimNCVTruth(scenario.numDimen,scenario.initState,scenario.qtrue,scenario.tsample,scenario.nsamples,scenario.observRegion);
        lastTrue(recycleMCdraws.runIdx).truesys    = truesys;
        lastTrue(recycleMCdraws.runIdx).H          = H;
        lastTrue(recycleMCdraws.runIdx).G          = G;
        lastTrue(recycleMCdraws.runIdx).qtrue      = scenario.qtrue;
        lastTrue(recycleMCdraws.runIdx).nsamples   = scenario.nsamples;
        lastTrue(recycleMCdraws.runIdx).tsample    = scenario.tsample;
        lastTrue(recycleMCdraws.runIdx).numDimen   = scenario.numDimen;
        lastTrue(recycleMCdraws.runIdx).numSensors = scenario.numSensors;
        lastTrue(recycleMCdraws.runIdx).initState  = scenario.initState;
        
    elseif ~isempty(recycleMCdraws) &&...
            recycleMCdraws.reuseNoise
        [truesys, H, G]     = deal(lastTrue(recycleMCdraws.runIdx).truesys, lastTrue(recycleMCdraws.runIdx).H, lastTrue(recycleMCdraws.runIdx).G);
    else
        [truesys, H, G]     = deal(lastTrue(1).truesys, lastTrue(1).H, lastTrue(1).G);
    end    
    
    save([scenario.simDataSaveDirectory, '/lastTrueFile'],'lastTrue')
elseif scenario.newTrajEachMCrun && (isempty(recycleMCdraws) || (~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise))
    [truesys, H, G]     = generateNdimNCVTruth(scenario.numDimen,scenario.initState,scenario.qtrue,scenario.tsample,scenario.nsamples,scenario.observRegion);
    
    lastTrueTemp.truesys    = truesys;
    lastTrueTemp.H          = H;
    lastTrueTemp.G          = G;
    lastTrueTemp.qtrue      = scenario.qtrue;
    lastTrueTemp.nsamples   = scenario.nsamples;
    lastTrueTemp.tsample    = scenario.tsample;
    lastTrueTemp.numDimen   = scenario.numDimen;
    lastTrueTemp.numSensors = scenario.numSensors;
    lastTrueTemp.initState  = scenario.initState;
    if ~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise
        lastTrue(recycleMCdraws.runIdx) = lastTrueTemp;
    else
        lastTrue = lastTrueTemp;
    end
    save([scenario.simDataSaveDirectory, '/lastTrueFile'],'lastTrue')
    
else
    [truesys, H, G]     = generateNdimNCVTruth(scenario.numDimen,scenario.initState,scenario.qtrue,scenario.tsample,scenario.nsamples,scenario.observRegion);
    
    lastTrueTemp.truesys    = truesys;
    lastTrueTemp.H          = H;
    lastTrueTemp.G          = G;
    lastTrueTemp.qtrue      = scenario.qtrue;
    lastTrueTemp.nsamples   = scenario.nsamples;
    lastTrueTemp.tsample    = scenario.tsample;
    lastTrueTemp.numDimen   = scenario.numDimen;
    lastTrueTemp.numSensors = scenario.numSensors;
    lastTrueTemp.initState  = scenario.initState;
    if ~isempty(recycleMCdraws) && recycleMCdraws.reuseNoise
        lastTrue(recycleMCdraws.runIdx) = lastTrueTemp;
    else
        lastTrue = lastTrueTemp;
    end
    save([scenario.simDataSaveDirectory, '/lastTrueFile'],'lastTrue')
end
trueSys = truesys;
trueSys.H = H;

if ~scenario.useLastTruth
    figure(13);
    plot(truesys.xtrue(1,:),truesys.xtrue(2,:),'-x','LineWidth',2)
    title('Possible truth trajectory okay (zero (0) to cancel)');
    axis(scenario.observRegion(1:4))

    if ~scenario.useFirstNew
    cont = input('enter 0 to cancel:  ');
    if cont ==0
        keyboard;
    end
    end
elseif scenario.plotExtraSetupInfo
    figure(13);
    plot(trueSys.xtrue(1,:),trueSys.xtrue(2,:),'-x','LineWidth',2)
    title('Possible truth trajectory okay (zero (0) to cancel)');
    axis(scenario.observRegion(1:4))

    pause(.5)
end
