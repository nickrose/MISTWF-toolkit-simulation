function [c,ceq] = constraints_eh_lin(p)
% constraints_eh_lin.m - this function it meant to be used in conjunction with
% the energy harvesting optimization simulations for constrains linearly
% related to the harvested energy, i.e. for the sum log objective
% 
% USAGE:
%    [c,ceq] = constraints_eh_lin(p)
%
% Where 'p' is the input vector of optimization variables and 'c' and 'ceq'
% are the constraint inequalities and equalities, respectively.
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

global sensorData

N = sensorData.numSensors;
d = sensorData.numElem;
H = sensorData.numEpochs;

% alpha2 = ( (repmat(sensorData.commParam.pathloss',1,H)./sensorData.commParam.h)*...
%             log(1./sensorData.commParam.PrBitError)...
%           ).^2;

if length(p) == N*d*H
    numVarMod = d;
else
    numVarMod = 1;
end
% expandedAlpha = reshape(repmat(alpha(:)',numVarMod,1),numVarMod*N,H);
% expandedAlpha2 = reshape(repmat(alpha2(:)',numVarMod,1),numVarMod*N,H);

powerUsed = reshape(p,N*numVarMod,H);
energyUsed = powerUsed .* repmat(sensorData.currentInterval,N*numVarMod,1);
%lossXdynLevels = expandedAlpha2 .* reshape (q,N*numVarMod,H);%expandedAlpha*...

% the total power constraint per time step
% if N > 1
    totPowConstr = sum(energyUsed,1) - sensorData.commParam.PtotConstr*ones(1,H);
% else
%     totPowConstr = [];
% end

% the summation of energy spent per time step
sumPL = sum(reshape(energyUsed,numVarMod,N*H),1);
sumPL = reshape(sumPL,N,H);

EminusPL =  sensorData.energyHarvestForOptimz - sumPL;

%combinedConstr = zeros(N*H,1);
% consecutive summation constraints; power can move forward from the
% current time step to the next (up to batt max); cannot use more energy
% than received; these are satisfied for enlarging summations
%
% the simplified constraint here accounts for causality and battery
% limitations
for jj = 1:H
    % simplified constraint
    causalityConstr(1+(jj-1)*N:N*jj) = sum(EminusPL(:,1:jj),2)';
    batteryConstr(1+(jj-1)*N:N*jj) = (sensorData.Emax - sum(EminusPL(:,1:jj),2))';
    sensorData.batteryStateOptmz(:,jj) = sum(EminusPL(:,1:jj),2);
end

c = [totPowConstr';...      sum_n p_nk < ptot
     -causalityConstr';...  sum_k=1^j p_nk T < sum_k=1^j E_nk   
     -batteryConstr';...     sum_k=1^j E_nk - sum_k=1^j p_nk < Emax 
     -p]; % positivity of dynamic number of levels

ceq = [];

if sensorData.plotOptmzProgress
    recordProgress(totPowConstr',[causalityConstr';batteryConstr'],p);
end

% =========================================================================
function recordProgress(totPowConstr,combinedConstr,q)
global progress

N = progress.numSensors;
H = progress.horizon;
numVarMod = length(q)/(N*H);


if isempty(progress.quantzLevels)
    progress.obj(1) =  abs(progress.objFunc(q));
    progress.quantzLevels(:,:,1) = reshape(sum(reshape(q,numVarMod,N*H),1),N,H);
    progress.constrValue(:,1) = [norm(totPowConstr);...
        min(abs(totPowConstr));...
        norm(combinedConstr);...
        min(abs(combinedConstr));...
        norm(-q);...
        ];
    progress.constrViol(:,1) = [1 * any(totPowConstr > 0);...
        2 * any(combinedConstr > 1e-8);...
        3 * any(-q > 1e-8)...
        ];
else
    progress.obj(end+1) =  progress.objFunc(q);
    progress.quantzLevels(:,:,end+1) = reshape(sum(reshape(q,numVarMod,N*H),1),N,H);
    progress.constrValue(:,end+1) = [norm([0;totPowConstr(totPowConstr > 0)]);...
        min(abs([0;totPowConstr(totPowConstr > 0)]));...
        norm([0;combinedConstr(combinedConstr > 1e-8)]);...
        min(abs([0;combinedConstr(combinedConstr > 1e-8)]));...
        norm(-q);...
        ];
    progress.constrViol(:,end+1) = [1 * any(totPowConstr > 0);...
        2 * any(combinedConstr > 1e-8);...
        3 * any(-q > 1e-8)...
        ];
end
if progress.onlinePlotting
    ColorOrder = progress.ColorOrder;

    numIter = length(progress.obj);
    % draw the objective and norm of the variables
    figure(11)
    semilogy(1:numIter, abs(progress.obj),'-k+');hold on;
    for nn = 1:N
        levels = reshape(progress.quantzLevels(nn,:,:),H,numIter);
        plot(1:(numIter), levels','Color',ColorOrder(1+mod(nn-1,21),:));hold on;
    end
    hold off;
    legend('MSE obj','Power Levels per node')
    title('Objective and variable values')
    
    % Constraint violations graphs
    figure(12)
    subplot(2,1,1)
    semilogy(repmat([1:numIter]',1,5), abs((progress.constrValue)'));
    legend('< Ptot','nrst < Ptot','Combined','nrst Combined','Pos semi-def var','location','west')
    title('Constraint satisfaction/violation')
    subplot(2,1,2)
    plot(repmat([1:numIter]',1,3), progress.constrViol');
    legend('< Ptot','Combined Constr','Pos semi-def var','location','west')
    xlabel('Viol: if value \geq 1 => at least one violation')
    thisAx = axis; thisAx(3:4) = [-1 5];
    axis(thisAx);
    % keyboard
    pause(.05)
end
