function metricsMC_eh(nameModifier,scenario,Metrics,aDate,here)
% metricsMC_eh.m - plots results for Monte Carlo runs of the optimization
% simulations for the energy harvesting scenarios.
%
% USAGE:
%   metricsMC_eh(nameModifier,scenario,Metrics,aDate,here)
%  
% INPUTS:
%   Metrics     
%
%   scenario      the scenario data; contains all the specifications
%                 sensing, communications, for how the data was generated,
%                 etc. 
%   nameModifier  a string to add to the filename when saving data and plots
%   aDate         (optional) an optional date to force into the filename,
%                 otherwise the current date is used
%   here          (optional) the folder to save the data and plots in, the
%                 default is the 'newResults' folder and requires the using
%                 to be in the main ./Code/ folder; can also simply set to
%                 TRUE to save in the current folder
%
% NO OUTPUTS (except for plots created)
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
% By: Nick Roseveare, Jan 2011

exitEarly=0;
if nargin > 3
    thedate = aDate;
elseif nargin == 0
    help metricsMC_eh
    exitEarly=1;
elseif isfield(scenario,'date')
    thedate = scenario.date;
else
    thedate = datestr(date,'yyyy-mm-dd');
end
if nargin > 4
    folderhere = here;
else
    folderhere = 0;
end

% if scenario.singleShot
%     if isempty(nameModifier)
%         nameModifier = 'singleShot';
%     else
%         nameModifier = [nameModifier '_singleShot'];
%     end
% end
LineStyleOrderWithMarker = {'-+','--^','-.d',':o','-*','--x', '-.s',':p'};
ColorOrder = ...
   [0    0     1;...
    0   .5     0;...
    1    0     0;...
    .5   0     1;...
    .75  0   .45;...
    .75 .75    0;...
    .25 .25  .25];

ColorOrder = [ColorOrder; ColorOrder*diag([.8 .7 .5]); ColorOrder*diag([.5  .4 .8]) ];


% numVarMod = scenario(1,1).optFormSetup.pMetricSize; % can be 1 or numDim
% H = scenario.optFormSetup.horizon;
% numSensors = scenario.numSensors;
% sizeMetrics = length(Metrics);
if ~exitEarly
% Save the data and possibly make plots
if ischar(folderhere)
    saveFolder = folderhere;
elseif folderhere
    saveFolder = './';
else
    saveFolder = '../newResults/';
end
if scenario.saveDataAndPlots
    if isfield(scenario,'fixNoiseSamples') && scenario.fixNoiseSamples
        reusedNoise = '_reuseNoise';%reusedNoise = 'reuseNoise';
    else
        reusedNoise = '';    
    end
    saveDataFilename = sprintf('%s/%srunData_%s_vs%s_%dMCruns%s',saveFolder,thedate,nameModifier,Metrics.varyOverStr,Metrics.numMC,reusedNoise);
    if scenario.checkForOW && exist(saveDataFilename,'file')
        if strcmp(input('The data file already exists, overwrite? (y/n)'),'y')
            save(saveDataFilename,'Metrics','scenario');
        end
    else
        save(saveDataFilename,'Metrics','scenario');
    end
end

% =========================== Make plots ===============================

if ~isfield(scenario,'makePlots') ||(isfield(scenario,'makePlots') && scenario.makePlots)

    h1 = figure(1);
    hold off;
    % MC runs index on rows, scenario variable indexed on columns
    sumRatePerf = [];
    validIdx = [];
    legendStr = {};
    methodList = fields(Metrics);
    for mtdIdx = 1:length(methodList)
        mtdStr = methodList{mtdIdx};
        if ~any(strcmp(mtdStr,{'scenario','numMC','varyOver','varyOverStr'}))
            validIdx(end+1) = mtdIdx;%#ok
            sumRatePerf(end+1,:) = inf*ones(1,length(Metrics.varyOver));%#ok
            if Metrics.(mtdStr)(1,1).sensorData.useOracle
                for vv = 1:length(Metrics.varyOver)
                    sumRatePerf(end,vv) = mean(abs([Metrics.(mtdStr)(:,vv).objectiveVal]));%#ok
                end
            else
                for vv = 1:length(Metrics.varyOver)
                    temp = [];
                    for mm = 1:Metrics.numMC
                        temp(mm) = sum(abs([Metrics.(mtdStr)(mm,vv).objectiveVal]));%#ok
                    end
                    sumRatePerf(end,vv) = mean(temp);%#ok
                end
            end
        end
        % Construct the legend
        switch mtdStr
            case 'stdOptimizRoutine'
                legendStr = ['UPBD' legendStr]; %#ok - put Upper bound first
            case 'predictedLayeredWF'
                legendStr{end+1} = 'MIST-WF-MPC'; %#ok     
            case 'noncausalknownLayeredWF'
                legendStr{end+1} = 'MIST-WF-NC'; %#ok
            case 'rfaRoutine'
                legendStr{end+1} = 'RFA'; %#ok
            otherwise
                if ~any(strcmp(mtdStr,{'scenario','numMC','varyOver','varyOverStr'}))
                    error('A method type was recorded in the Metrics struct which has not been defined')
                end
        end
        
    end
    plotIdx = 1;
    % Put uppper bound first
    reorderIdx = find(strcmp('stdOptimizRoutine',methodList(validIdx)));
    reorderIdx = [reorderIdx setdiff(1:length(validIdx),reorderIdx)];
    for mtdIdxLin = 1:length(methodList(validIdx))
        mtdIdx = reorderIdx(mtdIdxLin);
        if ~any(strcmp(methodList{validIdx(mtdIdx)},{'scenario','numMC','varyOver','varyOverStr'}))
            
%             ylabel('Sum Rate (nats)')
            ylabel('Inverse MSE','FontSize',14)
            switch Metrics.varyOverStr
                case 'nsamples'
                    if scenario.titleBars
                        title('Sum Rate vs. Number of Epochs','FontSize',14);
                    end
                    xlabel('No. of Epochs','FontSize',14)
                case 'numSensors'
                    if scenario.titleBars
                        title('Sum Rate vs. Number of Sensors','FontSize',14);
                    end
                    xlabel('No. of Sensors','FontSize',14)
                case 'dopplerFreqMax'
                    if scenario.titleBars
                        title('Sum Rate vs. Max Doppler Spread Time Step Product','FontSize',14);
                    end
                    xlabel('Max Doppler Spread Time Step Product','FontSize',14)
                    Metrics.varyOver = Metrics.varyOver .*[Metrics.scenario.tsample];
            end
            plot(Metrics.varyOver,sumRatePerf(reorderIdx(plotIdx),:),LineStyleOrderWithMarker{plotIdx},'MarkerSize',4*(length(legendStr)-plotIdx)+3,'LineWidth',.5*(length(legendStr)-plotIdx)+5,'Color',ColorOrder(plotIdx,:))
            hold on;
            plotIdx = plotIdx + 1;
         end
    end
    legend(legendStr,'Location','NorthWest','FontSize',14)
    hold off;
    
    if scenario.saveDataAndPlots
        saveName = sprintf('%s/%s_%s_vs%s_%dMCruns',saveFolder,thedate,nameModifier,Metrics.varyOverStr,Metrics.numMC);
        saveas(h1,saveName, 'fig');
        saveas(h1,saveName, 'psc2');
    end
    
    
%     h13=figure(13);hold off;
%     plot(scenario.trueSys.xtrue(1,:),scenario.trueSys.xtrue(2,:),'-o','LineWidth',2);
%     hold on;
%     plot(scenario.sensorPositions(1,:),scenario.sensorPositions(2,:),'rx','MarkerSize',10);hold on;
%     offset = scenario.observRegion([2 4])-scenario.observRegion([1 3]);
%     for ii = 1:numSensors
%         text(scenario.sensorPositions(1,ii)+offset(1)/40,scenario.sensorPositions(2,ii)+offset(2)/40,num2str(ii));hold on;
%     end
%     if scenario.titleBars
%         title('Sensor positions and Object trajectory');
%     end
%     axis(scenario.observRegion(1:4))
%     axis equal
%     hold off;
%     ylabel('North (du)')
%     xlabel('West (du)')
%     axis equal;
%     
%     if scenario.saveDataAndPlots
%         saveName = sprintf('%s/%s_%s_truthTrajectory',saveFolder,thedate,nameModifier);
%         saveas(h13,saveName, 'fig');
%         saveas(h13,saveName, 'psc2');
%     end
%     
%     sensorLgd = [];
%     for ii = 1:numSensors
%         sensorLgd{ii} = ['Node ' num2str(ii)];%#ok
%     end
%     
%     if scenario.singleShot
%         h1=figure(1);
%         for mm = 1:sizeMetrics
%             plot(Metrics(mm).bits','Color',ColorOrder(mm,:),'LineWidth',sizeMetrics-mm+1);hold on;
%         end
%         if scenario.titleBars
%             title('Bit usage')
%         end
%         legend(sensorLgd,'location','best');
%         axis auto
%         thisaxis = axis;
%         axis([thisaxis(1:2) -1e-3 thisaxis(4)])
%         xlabel('Time (tu)'),
%         ylabel('bits')
%         
%         if scenario.saveDataAndPlots
%             saveName = sprintf('%s/%s_%s_bitsUsage',saveFolder,thedate,nameModifier);
%             saveas(h1,saveName, 'fig');
%             saveas(h1,saveName, 'psc2');
%         end
%         hold off;
%         
%         
%         
%         h2=figure(2);
%         for mm = 1:sizeMetrics
%             plot(Metrics(mm).energyUsed','Color',ColorOrder(mm,:),'LineWidth',sizeMetrics-mm+1);hold on;
%         end
%         if scenario.titleBars
%             title('Energy usage')
%         end
%         axis auto
%         thisaxis = axis;
%         axis([thisaxis(1:2) 0 thisaxis(4)])
%         legend(sensorLgd,'location','best');
%         xlabel('Time (tu)'),
%         ylabel('Energy (eu)')
%         
%         if scenario.saveDataAndPlots
%             saveName = sprintf('%s/%s_%s_energyUsage',saveFolder,thedate,nameModifier);
%             saveas(h2,saveName, 'fig');
%             saveas(h2,saveName, 'psc2');
%         end
%         hold off;
%         
%         
%         
%         h3=figure(3);
%         for mm = 1:sizeMetrics
%             plot((Metrics(mm).sensorData.commParam.hHistoryOracle./repmat(Metrics(mm).sensorData.commParam.pathloss',1,H))','Color',ColorOrder(mm,:),'LineWidth',sizeMetrics-mm+1);hold on;
%         end
%         if scenario.titleBars
%             title('Actual Fading / Pathloss')
%         end
%         axis auto
%         thisaxis = axis;
%         axis([thisaxis(1:2) 0 thisaxis(4)])
%         legend(sensorLgd,'location','best');
%         xlabel('Time (tu)'),
%         ylabel('Fading/Pathloss')
%         
%         if scenario.saveDataAndPlots
%             saveName = sprintf('%s/%s_%s_fadingPathloss',saveFolder,thedate,nameModifier);
%             saveas(h3,saveName, 'fig');
%             saveas(h3,saveName, 'psc2');
%         end
%         hold off;
%         
%         
%         
%         h4=figure(4);
%         for mm = 1:sizeMetrics
%             plot(Metrics(mm).sensorData.eHistoryOracle','Color',ColorOrder(mm,:),'LineWidth',sizeMetrics-mm+1);hold on;
%         end
%         if scenario.titleBars
%             title('Energy harvested')
%         end
%         legend(sensorLgd,'location','best');
%         axis auto
%         thisaxis = axis;
%         axis([thisaxis(1:2) 0 thisaxis(4)])
%         xlabel('Time (tu)'),
%         ylabel('Energy arrival (eu)')
%         
%         if scenario.saveDataAndPlots
%             saveName = sprintf('%s/%s_%s_energyHarvested',saveFolder,thedate,nameModifier);
%             saveas(h4,saveName, 'fig');
%             saveas(h4,saveName, 'psc2');
%         end
%         hold off;
%         
%         
%         
%         h5=figure(5);
%         hold off;
%         for mm = 1:sizeMetrics
%             for nn = 1:numSensors
%                 stairs((0:H+1)*scenario.tsample,...
%                     [Metrics(mm).sensorData.energyRem(nn,:) Metrics(mm).batteryState(nn,:) Metrics(mm).batteryState(nn,end)],...
%                     LineStyleOrderWithMarker{1+mod(nn-1,8)},...
%                     'Color',ColorOrder(mod(nn+numSensors*(mm-1),7),:),'LineWidth',2*(sizeMetrics*numSensors+1-nn-numSensors*(mm-1))-1);%6*(numSensors+1-nn)/numSensors)
%                 hold on;
%             end
%         end
%         if scenario.titleBars
%             title('Battery state')
%         end
%         legend(sensorLgd,'location','best');
%         axis auto
%         thisaxis = axis;
%         axis([thisaxis(1:2) 0 thisaxis(4)])
%         xlabel('Time (tu)'),
%         ylabel('Battery level (eu)')
%         
%         if scenario.saveDataAndPlots
%             saveName = sprintf('%s/%s_%s_batteryLevels',saveFolder,thedate,nameModifier);
%             saveas(h5,saveName, 'fig');
%             saveas(h5,saveName, 'psc2');
%         end
%         hold off;
%         
%         
%         
%         h6=figure(6);
%         for mm = 1:sizeMetrics
%             plot(reshape(sum(reshape(Metrics(mm).sensorData.pHistoryOracle,numVarMod,numSensors*H),1),numSensors,H)','Color',ColorOrder(mm,:),'LineWidth',sizeMetrics-mm+1);hold on;
%         end
%         if scenario.titleBars
%             title('Local state uncertainty')
%         end
%         legend(sensorLgd,'location','best');
%         axis auto
%         thisaxis = axis;
%         axis([thisaxis(1:2) 0 thisaxis(4)])
%         xlabel('Time (tu)'),
%         ylabel('MSE')
%         
%         if scenario.saveDataAndPlots
%             saveName = sprintf('%s/%s_%s_localStateUnc',saveFolder,thedate,nameModifier);
%             saveas(h6,saveName, 'fig');
%             saveas(h6,saveName, 'psc2');
%         end
%         hold off;
%     else
%         1;
%     end    
% end
% 
% if scenario.singleShot && scenario.plotOptmzProgress && isfield(Metrics,'progress')
%     progress = Metrics.progress;
%     N = progress.numSensors;
%     H = progress.horizon;
%     
%     numIter = length(progress.obj);
%     % draw1the objective and norm of the variables
%     h11 = figure(11);
%     semilogy(1:numIter, abs(progress.obj),'-k+');hold on;
%     for nn = 1:N
%         levels = reshape(progress.quantzLevels(nn,:,:),H,numIter);
%         levels(levels > -1e-6 & levels < 0) = 0;
%         plot(1:(numIter), levels','Color',ColorOrder(1+mod(nn-1,21),:));hold on;
%     end
%     hold off;
%     legend('MSE obj','Quantz Levels per node')
%     title('Objective and variable values')
%     
%     if size(progress.constrValue,1) > 5
%         % Constraint violations graphs
%         h12 = figure(12);
%         subplot(2,1,1)
%         semilogy(repmat([1:numIter]',1,7), abs((progress.constrValue)'));
%         legend('< Ptot','nrst < Ptot','Causality','nrst Causality','Max Battery','nrst Max Battery','Pos semi-def var','location','west')
%         title('Constraint satisfaction/violation')
%         subplot(2,1,2)
%         plot(repmat([1:numIter]',1,4), progress.constrViol');
%         legend('< Ptot','Causality','Max Battery','Pos semi-def var','location','west')
%         xlabel('Viol: if value \geq 1 => at least one violation')
%         thisAx = axis; thisAx(3:4) = [-1 5];
%         axis(thisAx);
%     else
%         
%         % Constraint violations graphs: simplified constraint
%         h12 = figure(12);
%         subplot(2,1,1)
%         semilogy(repmat([1:numIter]',1,5), abs((progress.constrValue)'));
%         legend('< Ptot','nrst < Ptot','Combined','nrst Combined','Pos semi-def var','location','west')
%         title('Constraint satisfaction/violation')
%         subplot(2,1,2)
%         plot(repmat([1:numIter]',1,3), progress.constrViol');
%         legend('< Ptot','Combined Constr','Pos semi-def var','location','west')
%         xlabel('Viol: if value \geq 1 => at least one violation')
%         thisAx = axis; thisAx(3:4) = [-1 5];
%         axis(thisAx);
%         % keyboard
%         pause(.05)
%     end
%     if scenario.saveOptmzProgress && scenario.saveDataAndPlots
%         saveName = sprintf('%s/%s_%s_Optimiz-ObjAndVar',saveFolder,thedate,nameModifier);
%         saveas(h11,saveName, 'fig');
%         saveas(h11,saveName, 'psc2');
%         
%         saveName = sprintf('%s/%s_%s_Optimiz-constraints',saveFolder,thedate,nameModifier);
%         saveas(h12,saveName, 'fig');
%         saveas(h12,saveName, 'psc2');
%     end
% end

end
end