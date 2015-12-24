function [r, varargout] = generalMISTWF(margUtil,margUtilInv,params,energy,battmax,PtotMax,misc)
% layeredDirectionalSTWF.m - execute the layered directional water-filling
% algorithm for the given energy and initial marginal utility (zero 
% allocation mariginal utility).
%
% USAGE:
%       [r, exitflag] = generalMISTWF(margUtil,margUtilInv,params,energy,battmax,PtotMax,misc)
%
% DESCRIPTION:
% This algorithm solves the following problem (optimally/approx?)
%
%    max. sum_n sum_k a_nk * U(b_nk*r_nk+c_nk)
%     subject to
%       sum_n   d_nk < PtotMax
%       sum_k^j r_nk d_nk <= sum_k^j E_nk  forall j,n  
%       sum_k^j E_nk - sum_k^j r_nk d_nk <= battmax   forall j,n,
%       sum_n r_nk d_k <= P_tot  forall n.
%
% The distinction of this algorithm from the original layered WF algorithm
% is the use of 2-D data. This represents water filling across sensors and
% in time (n and k, respectively).
%
% The matrix quantities here use row indices for sensors, and column
% indices for time, i.e., input matrices will be NxK (N: no. sensors, K:
% no. of time steps).
%
% INPUT:
%   margUtil    the marginal utility function of the indiv. utility, i.e.,
%               the overall utility is defined in terms of this function as
%                  U_T = sum  a_nk*U(b_nk*r_nk + c_nk)  s.t.   sum d_nk*r_nk <= Q_T
%               and the marginal utility is V = dU/dr
%   margUtilInv the inverse of the marginal utility function
%   energy (E)  the amount of energy available from the harvester at each
%               time instance, a matrix quantity (or the avail. resources)
%   params      contains subfields for the scaling parameters a_nk, b_nk,
%               c_nk, d_nk which are each N x K.
%   battmax     the maximum battery capacity
%
%   misc        a dummy field for including various misc control variables
%
% OUPUT:
%   r           the allocated powers subject to the directionality
%               constraints
%
%   exitflag    (optional) informs the user when the algorithm exited
%               without an excessive number of "check" iterations
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
% By: Nick Roseveare, February 2012

% Defaults
playBack = 0;
plotProgress = 0;
plotExtraPriceResourceInfo = 0;
iterationsMax = 10;
elimNodeEpochThres = 1e-8;
% Optional override default settings
if nargin > 5
    strs = fields(misc);
    for ss = 1:length(strs);
        eval(sprintf('%s = misc.%s;',strs{ss},strs{ss}));
    end
end


[N,K] = size(params.a);
%L = intervals;%ones(1,n); % default value until algorithm needs to include the effects of the epoch
r = zeros(N,K);

if exist('W','var') && min(size(W))== 1 % use for the marginal utility in distrb estimation problem
    % Takes n*d x 1 and makes it n*d x 
    W = repmat(W,1,K);
    %W = W(:);
    %Wlin = W(:);
elseif ~exist('W','var')
    W = zeros(N,K);
end

% Initialize th algorithm parameters
k = 1;
T{k} = [];% in this two dimensional case elements of Tk become 2-tuples
V{k} = 1:N*K;
[VtimeStep,VsensorNode] = meshgrid(1:K,1:N);
VtimeStepIdx = VtimeStep(:);
VsensorNodeIdx = VsensorNode(:);
[time,node] = meshgrid(1:N:N*K,0:N-1);
VlinearIdx = time + node;
Etotal = zeros(N,K);
energyToUse = energy;

for ii = 1:K
    Etotal(:,ii) = sum(energy(:,1:ii),2);%sum(energy(:,1:ii)./repmat(L(1:ii),N,1),2);%  %****%
end
remainingPowerTotal = ones(K,1)*PtotMax;
AllTkSetSizeNodes = zeros(K,1);

blue = [0 0 1];
red  = [1 0 0];

% for investigating the shape of the resource and price functions
   % plotResAndPrices(r,margUtil,margUtilInv,params,W,Vplot,VsensorNodeIdx,VtimeStepIdx)

% ** METRIC ** plot initial loss profile
if plotProgress
figure(16)
maxlevel = max(max(energy))*1.05;
initMargUtil = priceFunction(r,margUtil,params,W,'all');
maxlevelLossStrt = max(max(1./initMargUtil))*1.05;
for nn = 1:N
    subplot(N,2,2*nn-1)
    bar(energy(nn,:),'stack')%,'FaceColor',blue)%,'bar_color','c')
    theAxis = axis;
    axis([theAxis(1:3) maxlevel]);
    if nn == 1
        title('Initial Available Energy')
    end
    ylabel(sprintf('Sensor %d',nn))
    if nn == N
        xlabel('Time steps')
    end
    figure(16)
    subplot(N,2,2*nn)
    bar(1./initMargUtil(nn,:),'stack')%,'FaceColor',red)
    theAxis = axis;
    axis([theAxis(1:3) maxlevelLossStrt]);
    if nn == 1
        title('Initial Marginal Utility profile')
    end
    if nn == N
        xlabel('Time steps')
    end
end
end
plotFlowPriceValues.margUtilIncr = inf;
plotFlowPriceValues.battPrice = ones(N,1)*inf;
plotFlowPriceValues.resAvailPrice = ones(N,1)*inf;
plotFlowPriceValues.epochTotPrice = ones(N,1)*inf;
plotFlowPriceValues.step = 0;

%************
% tStart = tic;
% Run the algorithm
Elast = Etotal;
iterations = 0;
stepCount = 0;
while any( any(Etotal > elimNodeEpochThres,1) & (remainingPowerTotal' > elimNodeEpochThres) ) && iterations < iterationsMax
    stepCount = stepCount + 1;
    garbageCollectSubsetIdx = [];
    
    %++++++++++++++++++++++++++++++++++++++++
    % First stage increment:
    %++++++++++++++++++++++++++++++++++++++++
    % ========= Find minimum increment from marginal utility ==========
    for k = 1:length(T)
        % Find new difference levels
        marginal = priceFunction(r,margUtil,params,W,'all');
        marginal = marginal(:);
        m = max(marginal(V{k}));
        T{k} = V{k}(abs(m - marginal(V{k}))<elimNodeEpochThres);
        nextLowestPrice = max(marginal(setdiff(V{k},T{k})));
        if isempty(nextLowestPrice)
            d(k) = -1; %#ok
        else
            d(k) = nextLowestPrice;%#ok
        end
    end
    % ================================================================= 
	
    % For debugging
    if plotExtraPriceResourceInfo
    plotResAndPrices([],[],margUtil,margUtilInv,params,W,[],VsensorNodeIdx,VtimeStepIdx);
    end
    
    %++++++++++++++++++++++++++++++++++++++++
    % Second stage constraint increments:
    %++++++++++++++++++++++++++++++++++++++++
    % ===== Find total-epoch constraint remaining energy flow =====
    commonTotFlowPrice = zeros(1,K);
    %    tStepsHere = unique(VtimeStepIdx([T{:}]))';
    TnonUnique = [T{:}];
    allSetIdx = unique(TnonUnique);
    NodesThisTimeStep = {};
    AllTkSetSizeNodes = zeros(1,K);
    for mm = allSetIdx
%         if length(NodesThisTimeStep) <= VtimeStepIdx(mm) || isempty(NodesThisTimeStep{VtimeStepIdx(mm)})
            NodesThisTimeStep{VtimeStepIdx(mm)} = find(VtimeStepIdx(TnonUnique) == VtimeStepIdx(mm))';%#ok
            AllTkSetSizeNodes(VtimeStepIdx(mm)) = length(NodesThisTimeStep{VtimeStepIdx(mm)}); % number of updated nodes in each time step
%         end
        RnEven = remainingPowerTotal(VtimeStep(mm))/AllTkSetSizeNodes(VtimeStepIdx(mm));
        for kk = NodesThisTimeStep{VtimeStepIdx(mm)}
            jj = TnonUnique(kk);
            commonTotFlowPrice(VtimeStepIdx(mm)) =...
                max(commonTotFlowPrice(VtimeStepIdx(mm)),...
                    priceFunction(r(VsensorNodeIdx(jj),VtimeStepIdx(jj))+RnEven,...
                    margUtil,params,W,VsensorNodeIdx(jj),VtimeStepIdx(jj))...
                   );
        end
        if plotExtraPriceResourceInfo
        plotResAndPrices(repmat(RnEven,N,K),[],margUtil,margUtilInv,params,W,TnonUnique(NodesThisTimeStep{VtimeStepIdx(mm)}),VsensorNodeIdx,VtimeStepIdx,1)%#ok
        end
    end

    % =====================================================================
    

    for k = 1:length(T)
        % Set next powers of bins/epochs
        if ~isempty(T{k}) %> 0

            TstepsThisNode = {};
            NumTstepsThisNode = zeros(N,1);
            for mm = T{k}
                TstepsThisNode{VsensorNodeIdx(mm)} = T{k}(find(VsensorNodeIdx(T{k}) == VsensorNodeIdx(mm))');%#ok
                IdxStepsThisNode{VsensorNodeIdx(mm)} = find(VsensorNodeIdx(T{k}) == VsensorNodeIdx(mm))';%#ok
                NumTstepsThisNode(VsensorNodeIdx(mm)) = length(TstepsThisNode{VsensorNodeIdx(mm)});
            end
                        
            commonEflowPrice = zeros(N,K);
            % ========== Find maximum energy availability ============
            for mm = T{k}
                % TthisNodeEflow{VsensorNodeIdx(mm)} = find(VsensorNodeIdx(T{k}) == VsensorNodeIdx(mm))';%#ok

%                 TstepsThisNode{VsensorNodeIdx(mm)} = find(VsensorNodeIdx(T{k}) == VsensorNodeIdx(mm))';%#ok
%                 NumTstepsThisNode(VtimeStepIdx(mm)) = length(TstepsThisNode{VsensorNodeIdx(mm)}); %#ok number of updated nodes in each time step
                
                RnEven = Etotal(VsensorNodeIdx(mm),VtimeStepIdx(mm))...
                                    /                    NumTstepsThisNode(VsensorNodeIdx(mm)); %
                                                           % length(TthisNodeEflow(find(TthisNodeEflow{VsensorNodeIdx(mm)}==mm):end))
                                                           %length(1:find(TstepsThisNode{VsensorNodeIdx(mm)} == mm));
                % Choose base price for all allocation functions - this
                % ensures feasibility, but will not acheive a active
                % constraint in a single step as before in the more simple
                % concave utility maximization
                    % for nn = TstepsThisNode{VsensorNodeIdx(mm)}
                        % jj = T{k}(nn);

                    commonEflowPrice(VsensorNodeIdx(mm),VtimeStepIdx(mm)) =...
                                    max([commonEflowPrice(VsensorNodeIdx(mm),VtimeStepIdx(mm)),...
                                        priceFunction(r(VsensorNodeIdx(mm),VtimeStepIdx(mm))+RnEven,...
                                        margUtil,params,W,VsensorNodeIdx(mm),VtimeStepIdx(mm))...
                                       ]); 
                    % end
                % For debugging:
                if plotExtraPriceResourceInfo
                plotResAndPrices(repmat(RnEven,N,K),[],margUtil,margUtilInv,params,W,T{k}(IdxStepsThisNode{VsensorNodeIdx(mm)}),VsensorNodeIdx,VtimeStepIdx,2)
                end

            end
            % =========================================================
            
            % ** METRIC **
            if plotProgress
            figure(17)
            % \/ only useful in removing small values for plotting when the
            % 'flow' is the actual resource, not the price
            %  commonEflowPrice(abs(commonEflowPrice)<elimNodeEpochThres) =
            %  0;
            maxlevel = max(max(commonEflowPrice))*1.05; 
            for nn = 1:N
                subplot(N+1,4,4*nn)
                
                bar(commonEflowPrice(nn,:)')%,'FaceColor',blue)
                theAxis = axis;
                axis([theAxis(1:3) maxlevel]);      
                if nn == 1
                    title('Causal Energy Available Price')
                end
                if nn == N
                    xlabel('Time steps')
                end
            end
            end
            
            % ===== Find availability of energy under storage constr =====
            batteryDebt = zeros(N,1);
            for mm = T{k}

                %TthisNode = find(VsensorNodeIdx(T{k}) == VsensorNodeIdx(mm));
                %TkSetSizeTime(VsensorNodeIdx(mm)) = length(TthisNode);% number of updated time step for each node
                %fillIndicesSoFarN = T{k}(TthisNodeEflow{VsensorNodeIdx(mm)}(1:find(T{k}(TthisNodeEflow{VsensorNodeIdx(mm)}) == mm)));

                % If some energy has already been drawn from another epoch,
                % we must account for this by a reduction in the available
                % (predicted) battery supply
                    % sensorNo = VsensorNodeIdx(mm);%VsensorNodeIdx(fillIndicesSoFarN);%unique(VsensorNodeIdx(fillIndicesSoFarN)); % <= this should be a single scalar, if not, then we have issues
                    % for nn = sensorNo'
                for nn = TstepsThisNode{VsensorNodeIdx(mm)}
                    idxToHere = find(TstepsThisNode{VsensorNodeIdx(nn)}==nn);
                    stepsUpToHere = TstepsThisNode{VsensorNodeIdx(mm)}(1:idxToHere);
                    batteryDebt(VsensorNodeIdx(nn)) = ...
                        min(batteryDebt(VsensorNodeIdx(nn)),...
                        sum(energy(VsensorNodeIdx(nn),VtimeStepIdx(stepsUpToHere)) ...
                        - params.d(VsensorNodeIdx(nn),VtimeStepIdx(stepsUpToHere)) ...
                        .* r(VsensorNodeIdx(nn),VtimeStepIdx(stepsUpToHere))...
                        )...
                        );
                 end
            end
            batteryDebtFlowPrice = zeros(N,K);
            for mm = T{k}
                if VtimeStepIdx(mm) > 1
                    RnEven = max((battmax + batteryDebt(VsensorNodeIdx(mm))),0) / NumTstepsThisNode(VsensorNodeIdx(mm));
                    batteryDebtFlowPrice(VsensorNodeIdx(mm),VtimeStepIdx(mm)) = ...
                        max([batteryDebtFlowPrice(VsensorNodeIdx(mm),VtimeStepIdx(mm)),...
                             priceFunction(r(VsensorNodeIdx(mm),VtimeStepIdx(mm))+RnEven,...
                             margUtil,params,W,VsensorNodeIdx(mm),VtimeStepIdx(mm))...
                            ]);
                end
                if plotExtraPriceResourceInfo
                plotResAndPrices(repmat(RnEven,N,K),[],margUtil,margUtilInv,params,W,T{k}(IdxStepsThisNode{VsensorNodeIdx(mm)}),VsensorNodeIdx,VtimeStepIdx,3);%#ok
                end
            end
            
            % ===========================================================
            
            % Debugging: plots prices on marginal distributions
            if plotExtraPriceResourceInfo
            plotResAndPrices([],repmat(commonTotFlowPrice,N,1),margUtil,margUtilInv,params,W,T{k},VsensorNodeIdx,VtimeStepIdx,1);
            plotResAndPrices([],repmat(max(commonEflowPrice,[],2),1,K),margUtil,margUtilInv,params,W,T{k},VsensorNodeIdx,VtimeStepIdx,2);
            plotResAndPrices([],repmat(max(batteryDebtFlowPrice,[],2),1,K),margUtil,margUtilInv,params,W,T{k},VsensorNodeIdx,VtimeStepIdx,3);
                pause(.3)
            end
            % =========== Find the minimum allocation increment =========

            % For debug - record price changes
             %   plotFlowPriceValues = plotFlowPrices(plotFlowPriceValues,stepCount,k,d(k),batteryDebtFlowPrice,commonEflowPrice,commonTotFlowPrice);
            
            rnew = zeros(N,K);
            for mm = T{k}
                if d(k) < 0
                    if VtimeStepIdx(mm) > 1
                        maxPrice = max([...
                                    batteryDebtFlowPrice(VsensorNodeIdx(mm),:),...
                                    commonEflowPrice(VsensorNodeIdx(mm),:),...
                                    commonTotFlowPrice(VtimeStepIdx(mm))]);
                    else
                        maxPrice = max([...
                                    commonEflowPrice(VsensorNodeIdx(mm),:),...
                                    commonTotFlowPrice(VtimeStepIdx(mm))]);      
                    end
                else
                    if VtimeStepIdx(mm) > 1
                        maxPrice = max([...
                                    d(k),...
                                    batteryDebtFlowPrice(VsensorNodeIdx(mm),:),...
                                    commonEflowPrice(VsensorNodeIdx(mm),:),...
                                    commonTotFlowPrice(VtimeStepIdx(mm))]);
                    else
                        maxPrice = max([...
                                    d(k),...
                                    commonEflowPrice(VsensorNodeIdx(mm),:),...
                                    commonTotFlowPrice(VtimeStepIdx(mm))]);
                    end
                end              
                rnew(VsensorNodeIdx(mm),VtimeStepIdx(mm)) = ...
                        max(r(VsensorNodeIdx(mm),VtimeStepIdx(mm)),...
                            allocResFunction(maxPrice,margUtilInv,params,W,VsensorNodeIdx(mm),VtimeStepIdx(mm))...
                           );
                
            end
            rinc = zeros(N,K);
            for mm = fliplr(T{k})
                rinc(VsensorNodeIdx(mm),VtimeStepIdx(mm)) = rnew(VsensorNodeIdx(mm),VtimeStepIdx(mm)) - r(VsensorNodeIdx(mm),VtimeStepIdx(mm));
                
                r(VsensorNodeIdx(mm),VtimeStepIdx(mm)) = rnew(VsensorNodeIdx(mm),VtimeStepIdx(mm));

                updateIdx = intersect(V{k}(VsensorNodeIdx(V{k})==VsensorNodeIdx(mm)), VlinearIdx(VsensorNodeIdx(mm),VtimeStepIdx(mm):K));
                Etotal(updateIdx) = Etotal(updateIdx) ...
                                    - params.d(VsensorNodeIdx(mm),VtimeStep(mm))*rinc(VsensorNodeIdx(mm),VtimeStepIdx(mm));
                
                energyToUse(VsensorNodeIdx(mm),VtimeStep(mm)) = ...
                                    energyToUse(VsensorNodeIdx(mm),VtimeStep(mm)) ...
                                        - params.d(VsensorNodeIdx(mm),VtimeStep(mm))*rinc(VsensorNodeIdx(mm),VtimeStep(mm));
            end 
            % =====================================================
        end
    end
    
        
    if ~isempty([V{:}])
        for mm = unique([T{:}])
            remainingPowerTotal(VtimeStep(mm)) = remainingPowerTotal(VtimeStep(mm)) ...
                                                 - params.d(VsensorNodeIdx(mm),VtimeStep(mm))*rinc(VsensorNodeIdx(mm),VtimeStep(mm));
        end
    end
    
    for k = 1:length(T)
        % Set next powers of bins/epochs
        if ~isempty(T{k}) > 0
            energyUsageProfile = energy - r.*params.d;
            splitIdxBorrowList = [];
            
            for nn = 1:N
                for ii = K:-1:2
                    if energyToUse(nn,ii) < 0
                        Etotal(nn,ii-1) = Etotal(nn,ii-1) + energyToUse(nn,ii);
                        energyToUse(nn,ii-1) = energyToUse(nn,ii-1) + energyToUse(nn,ii);
                        energyToUse(nn,ii) = 0;
                    end
                    if energyUsageProfile(nn,ii) < - elimNodeEpochThres
                        if abs(energyUsageProfile(nn,ii)) > battmax - elimNodeEpochThres
                            splitIdxBorrowList(end+1) = VlinearIdx(nn,ii); %#ok
                            % THIS IS A CRUX
                            % Not sure how to decompose problem since it
                            % is 2-D. How does one separate the beginning
                            % and end of the layering in time without
                            % separating the dependence at each time step
                            % on the total power constraint:
                            %  -> one Tk subset could possible not include
                            %  or account for the nodes from other Tk
                            %  subsets in the same time step?
                            % -> sol'n: create and update the total active
                            % number of node in each time step OUTSIDE of
                            % the loops with iterate over the Tk sets

                        end
                        energyUsageProfile(nn,ii-1) = energyUsageProfile(nn,ii-1) + energyUsageProfile(nn,ii);
                        energyUsageProfile(nn,ii) = 0;
                    end
                end
            end
            if any(sum(energy-params.d.*r,2) < - elimNodeEpochThres)
                 warning('Usage profile incorrect: more energy used than available')
            end

            
            % ** METRIC **
            if plotProgress
            figure(17)
            maxlevel = max(max(energy-r.*params.d))*1.05;             
            theMin = min(min(energy-r.*params.d));
            minlevel = sign(theMin)*(1+sign(theMin)*.2)*abs(theMin);             
            for nn = 1:N
                subplot(N+1,4,4*nn-1)
                energyUsePlot = energy(nn,:)-r(nn,:).*params.d(nn,:);
                energyUsePlot(abs(energyUsePlot) < elimNodeEpochThres) = 0;
                bar(energyUsePlot')%,'FaceColor',blue)
                theAxis = axis;
                axis([theAxis(1:2) minlevel maxlevel]);      
                if nn == 1
                    title('Energy Usage Difference')
                end
                
                if nn == N
                    xlabel('Time steps')
                end
            end
            end
            %************
            
            if all(all(Etotal == Elast))
                %keyboard
                iterations = iterations+1;
            else
                iterations = 0;
            end
            Elast = Etotal;
            
            %*** For Debugging ***
%             if toc(tStart) > 20
%                 keyboard
%             end
%             tStart = tic;
            
            % Eliminate epochs/bins whose remaining energy is zero
            elimIdx = find(Etotal <= elimNodeEpochThres);
            V{k} = setdiff(V{k},elimIdx);

            
            % Check for reaching maximum energy usage; indiv maximum usage,
            % and summed borrowing usage
            if ~isempty(V{k})
            for nn = 1:N
                thisNodeVidx = VlinearIdx(VsensorNode == nn);
                splitNodesThisSensor = intersect(splitIdxBorrowList,thisNodeVidx(:));
                splitIdxList = sort(splitNodesThisSensor,'ascend');%sort(unique([find(r(nn,:) >= battmax) splitNodesThisSensor]),'ascend');
                %                 r(r > battmax) = battmax;
                if ~isempty(splitIdxList);
                    knext = length(T) + 1;
                    splitIdxList = [1 splitIdxList];%#ok
                    EtotalPrev = 0;
                    for jj = 2:length(splitIdxList)
                        j = splitIdxList(jj);
                        jMinus1 = splitIdxList(jj-1);
                        % Split T into multiple subproblems, one on either side of
                        % the eliminated epochs/bins and one between each
                        V{knext} = intersect(thisNodeVidx(:)', V{k}(V{k} < j & V{k} >= jMinus1));
                        T{knext} = [];


                        % Remove the split problem sets from the current V
                        % and adjust energy resources availability list
                        if ~isempty(V{knext})
                            EtotalSubset = V{knext}(VsensorNodeIdx(V{knext})==VsensorNodeIdx(nn));
                            EtotalPrev = Etotal(EtotalSubset(end));
                            for rr = V{knext}
                                V{k}(V{k}==rr) = [];
                            end
                            % Update E total...resource in beginning epochs can
                            % no longer pass to later epochs
                            Etotal(V{k}(VsensorNodeIdx(V{k})==VsensorNodeIdx(mm))) = ...
                                Etotal(V{k}(VsensorNodeIdx(V{k})==VsensorNodeIdx(mm)))...
                                - EtotalPrev;
                        end


                        
                        knext = knext + 1;
                    end
                    % The original V{k} now contains all the nodes except
                    % for those time step for the nodes with battery
                    % maximal use acheived, i.e. something like
                    %
                    % updated original V:
                    %       1 4 7 10 13
                    %                14
                    %           9 12 15
                    
                end
            end
            end
        else
            garbageCollectSubsetIdx(end+1) = k; %#ok
        end
    end
    for jj = fliplr(garbageCollectSubsetIdx)
        T(jj) = [];
        V(jj) = [];
        %d(jj) = [];
    end
    
    % Eliminate epochs whose total power consumption has met the
    % per epoch power constraint
    for jj = K:-1:1
        if remainingPowerTotal(jj) <= elimNodeEpochThres
            for k = 1:length(V)
                % Set next powers of bins/epochs
                if ~isempty(V{k}) > 0
                    V{k}(VtimeStepIdx(V{k})==jj) = [];
                end
            end
        end
    end
    
    
    energyToUsePlot = energyToUse;
%     energyToUsePlot(abs(energyToUsePlot) < elimNodeEpochThres) = 0;
    
    % ** METRIC **
    if plotProgress
    figure(17)
     
    for nn = 1:N
        newMarg = priceFunction(r,margUtil,params,W,'all');
        maxlevelLoss = max(max(1./newMarg))*1.05; 
        subplot(N+1,4,4*nn-3)
        bar([1./(initMargUtil(nn,:)'), 1./newMarg(nn,:)'-1./initMargUtil(nn,:)'],'stack');%,'FaceColor',blue,'FaceColor',red)
        theAxis = axis;
        axis([theAxis(1:3) maxlevelLoss]);              
        if nn == 1
            title('Initial and Updated Marginal Utility Profile')
        end
        ylabel(sprintf('Sensor %d',nn))
        
        maxlevel = max(max(energyToUsePlot))*1.05; 
        subplot(N+1,4,4*nn-2)
        bar(energyToUsePlot(nn,:)')%,'FaceColor',blue)
        theAxis = axis;
        axis([theAxis(1:3) maxlevel]);        
        if nn == 1
            title('Energy Remaining')
        end
        if nn == N
            xlabel('Time steps')
        end
    end
    subplot(N+1,4,4*N+1)
    bar(remainingPowerTotal')%,'FaceColor',blue)
    ylabel('Total Constraint Slack')
    xlabel('Time steps')
    end
    %************
    
    if playBack
        pause%(.7)
    end
end
if nargout > 1
    if iterations < iterationsMax
        varargout{1} = 1;
    else
        varargout{1} = 0;
    end
end

% =========================================================================
function lambda = priceFunction(r,margUtil,params,W,nIdx,kIdx)
if ischar(nIdx) % if nIdx is 'all'
    lambda = params.a .*params.b ./params.d ...
        .* margUtil(params.b.*r+params.c,W);
else
    lambda = params.a(nIdx,kIdx).*params.b(nIdx,kIdx)./params.d(nIdx,kIdx) ...
        .* margUtil(params.b(nIdx,kIdx) .* r + params.c(nIdx,kIdx),W(nIdx,kIdx));
end

% =========================================================================
function r = allocResFunction(lambda,margUtilInv,params,W,nIdx,kIdx)

if ischar(nIdx) % if nIdx is 'all'
    r = 1./params.b .* ...
        (...max(
        margUtilInv(params.d ./(params.a .* params.b) .* lambda, W)...<- use index here to retrieve correct dynamic range (in case of distb estm problem)
            - params.c...
            );%,zeros(size(lambda)));
else
    r = 1./params.b(nIdx,kIdx) .* (...max(... 
                                      margUtilInv(params.d(nIdx,kIdx) ...
                                            ./(params.a(nIdx,kIdx) .* params.b(nIdx,kIdx)) ...
                                                  .* lambda, W(nIdx,kIdx)...<- use index here to retrieve correct dynamic range (in case of distb estm problem)
                                                  )...
                                     - params.c(nIdx,kIdx)...
                                   );%,0);
end
% % % For the case when we need to worry about keeping the res allocations positve %
% if ischar(nIdx) % if nIdx is 'all'
%     r = 1./params.b .* ...
%         (max(
%         margUtilInv(params.d ./(params.a .* params.b) .* lambda, W)...<- use index here to retrieve correct dynamic range (in case of distb estm problem)
%             - params.c...
%             ),zeros(size(lambda)));
% else
%     r = 1/params.b(nIdx,kIdx) .* (max(... 
%                                       margUtilInv(params.d ...
%                                             ./(params.a(nIdx,kIdx) .* params.b(nIdx,kIdx)) ...
%                                                   .* lambda, W(nIdx,kIdx)...<- use index here to retrieve correct dynamic range (in case of distb estm problem)
%                                                   )...
%                                      - params.c(nIdx,kIdx)...
%                                    ,0));
% end

% =========================================================================
function plotResAndPrices(r,price,margUtil,margUtilInv,params,W,Vplot,VsensorNodeIdx,VtimeStepIdx,markIdx)
    persistent minlevell maxlevell minlevelr maxlevelr figRes figPrice
    
%     dcm_obj = datacursormode(figRes);
%     set(dcm_obj,'UpdateFcn',@plotPriceInfoUpdateFunc)
    
    MarkerType = {'b+','r*','kd','go','cx', 'ms','bp'};
    if nargin <= 9
        markIdx = 1;
        replot = 1;
        figRes = figure(15);
        figPrice = figure(14);
%         dcm_objRes = datacursormode(figRes);
%         set(dcm_objRes,'UpdateFcn',@plotPriceInfoUpdateFunc)
%         dcm_objPrice = datacursormode(figPrice);
%         set(dcm_objPrice,'UpdateFcn',@plotPriceInfoUpdateFunc)
    else
        replot = 0;
    end
    res = 100;
    [N,K] = size(W);
    rLog =  logspace(-4,10,res);% logspace(log10(min([min(max(r,1e-20)), 1e-3]))-1,log10(max([max(r), 1e8]))+1,res);
    lambdaInputPrice = logspace(-25,1,res);%logspace(log10(min([min(max(price,1e-30)), 1e-25]))-1,log10(max([max(price), 1e3]))+1,res);
    lambdaOutputPlot = zeros(res,K,N);
    rAllocPlot = zeros(res,K,N);
    
if replot 
    warning('off','MATLAB:Axes:NegativeDataInLogAxis')
    for rr = 1:res
        rAll = repmat(rLog(rr),N,K);
        lamAll = repmat(lambdaInputPrice(rr),N,K);
        lambdaOutputPlot(rr,:,:) = priceFunction(rAll,margUtil,params,W,'all')';
        rAllocPlot(rr,:,:) = allocResFunction(lamAll,margUtilInv,params,W,'all')';
    end
    maxlevell = max(max(max(lambdaOutputPlot)))*1e2;
    maxlevelr = max(max(max(rAllocPlot)))*1e2;
    minlevell = min(min(min(lambdaOutputPlot)))*1e-2;
    minlevelr = max(min(min(min(rAllocPlot))),1e-1)*1e-2;

    for nn = 1:N
        figure(figRes)
        subplot(N,1,nn)
        loglog(repmat(rLog',1,K),lambdaOutputPlot(:,:,nn))
        grid on;
        ylabel(sprintf('Sensor %d',nn))
        axis([min(rLog) max(rLog) minlevell maxlevell]);
        if nn == 1
            title('V: price as function of resource')
        end
        figure(figPrice)
        subplot(N,1,nn)
        loglog(repmat(lambdaInputPrice',1,K),rAllocPlot(:,:,nn))
        grid on;
        ylabel(sprintf('Sensor %d',nn))
        axis([min(lambdaInputPrice) max(lambdaInputPrice) minlevelr maxlevelr]);
        if nn == 1
            title('V^{-1}: resource as function of price')
        end
        
    end
end
    if ~isempty(r)
        for ii = 1:length(Vplot)
            mm = Vplot(ii);
            figure(15);
            subplot(N,1,VsensorNodeIdx(mm))
            hold on;
            thePrice = priceFunction(r(mm),margUtil,params,W,VsensorNodeIdx(mm),VtimeStepIdx(mm));
            loglog(r(mm),thePrice,MarkerType{markIdx},'MarkerSize',8,'LineWidth',1.2)
            axis([min(rLog) max(rLog) minlevell maxlevell]);
            hold off;
        end
    end
    if ~isempty(price)
        for ii = 1:length(Vplot)
            mm = Vplot(ii);
            figure(14);
            subplot(N,1,VsensorNodeIdx(mm))
            hold on;
            resAll = allocResFunction(price(mm),margUtilInv,params,W,VsensorNodeIdx(mm),VtimeStepIdx(mm));
            loglog(price(mm),resAll,MarkerType{markIdx},'MarkerSize',8,'LineWidth',1.2)
            axis([min(lambdaInputPrice) max(lambdaInputPrice) minlevelr maxlevelr]);
            hold off;
        end
        
    end


% =========================================================================    
function plotFlowPriceValues = plotFlowPrices(plotFlowPriceValues,stepCount,...
                        k,d,batteryDebtFlowPrice,commonEflowPrice,commonTotFlowPrice)
persistent maxLevel minLevel
warning('off','MATLAB:Axes:NegativeDataInLogAxis')
if k == 1
if k > 1 && length(plotFlowPriceValues) < k
    plotFlowPriceValues(k).margUtilIncr  = d;
    plotFlowPriceValues(k).battPrice     = max(batteryDebtFlowPrice,[],2);
    plotFlowPriceValues(k).resAvailPrice = max(commonEflowPrice,[],2);
    plotFlowPriceValues(k).epochTotPrice = max(commonTotFlowPrice,[],2);
    plotFlowPriceValues(k).step          = stepCount;
else
    plotFlowPriceValues(k).margUtilIncr(end+1)    = d;
    plotFlowPriceValues(k).battPrice(:,end+1)     = max(batteryDebtFlowPrice,[],2);
    plotFlowPriceValues(k).resAvailPrice(:,end+1) = max(commonEflowPrice,[],2);
    plotFlowPriceValues(k).epochTotPrice(:,end+1) = max(commonTotFlowPrice,[],2);
    plotFlowPriceValues(k).step(end+1)            = stepCount;
end
N = size(plotFlowPriceValues(1).battPrice,1);

newPrices = [d(d>0); batteryDebtFlowPrice(:); commonEflowPrice(:); commonTotFlowPrice(:)];
newPrices = newPrices(~isinf(newPrices)&~isnan(newPrices));
if isempty(maxLevel)
    maxLevel = max(max(newPrices));
    minLevel = min(min(newPrices));    
else
    maxLevel = max([maxLevel, max(max(newPrices))]);
    minLevel = min([minLevel, min(min(newPrices))]);
end
thisAxis = [0 stepCount minLevel*(.9+.2*sum(minLevel < 0)) maxLevel*1.1];
figure(20+k)
subplot(4,1,1)
semilogy(plotFlowPriceValues(k).step,[plotFlowPriceValues(k).margUtilIncr])
axis(thisAxis);
title('Marginal Utility Increment Price')
subplot(4,1,2)
semilogy(repmat([plotFlowPriceValues(k).step]',1,N),[plotFlowPriceValues(k).battPrice]')
axis(thisAxis);
title('Battery Flow Price')
subplot(4,1,3)
semilogy(repmat([plotFlowPriceValues(k).step]',1,N),[plotFlowPriceValues(k).resAvailPrice]')
axis(thisAxis);
title('Available Resource Price')
subplot(4,1,4)
semilogy(repmat([plotFlowPriceValues(k).step]',1,N),[plotFlowPriceValues(k).epochTotPrice]')
axis(thisAxis);
title('Epoch Total Constraint Price')
pause(.5)
end
