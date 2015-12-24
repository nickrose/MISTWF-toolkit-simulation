function [c,ceq] = constraints_eh_simple(q)
% constraints_eh_simple.m - this function is meant to be used in conjunction with
% the energy harvesting optimization simulations.
% 
% USAGE:
%    [c,ceq] = constraints_eh_simple(q)
%
% Where 'q' is the input vector of optimization variables and 'c' and 'ceq'
% are the constraint inequalities and equalities, respectively.
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
% By: Nick Roseveare, December 2011

global sensorData

N = sensorData.numSensors;
d = sensorData.numElem;
H = sensorData.horizon;

alpha2 = ( (repmat(sensorData.commParam.pathloss',1,H)./sensorData.commParam.h)*...
            log(1./sensorData.commParam.PrBitError)...
          ).^2;
% alpha = ( (repmat(sensorData.commParam.pathloss',1,H)./sensorData.commParam.h)*...
%             log(1./sensorData.commParam.PrBitError)...
%           );

if length(q) == N*d*H
    numVarMod = d;
else
    numVarMod = 1;
end
% expandedAlpha = reshape(repmat(alpha(:)',numVarMod,1),numVarMod*N,H);
expandedAlpha2 = reshape(repmat(alpha2(:)',numVarMod,1),numVarMod*N,H);

powerUsed = expandedAlpha2 .* reshape (q,N*numVarMod,H);
lossXdynLevels = expandedAlpha2 .* reshape (q,N*numVarMod,H);%expandedAlpha*...

% the total power constraint per time step
totPowConstr = [sum(powerUsed,1) - sensorData.commParam.PtotConstr^2*ones(1,H)]';

% the summation of energy spent per time step
energySpentPerSensorPerTimeStep = sum(reshape(lossXdynLevels, numVarMod,N*H),1);
energySpentPerSensorPerTimeStep = reshape(energySpentPerSensorPerTimeStep,N,H);

EngyRcvdMinusSpent =  sensorData.energyHarvestForOptimz - energySpentPerSensorPerTimeStep;

combinedConstr = zeros(N*H,1);
% consecutive summation constraints; power can move forward from the
% current time step to the next (up to batt max); cannot use more energy
% than received; these are satisfied for enlarging summations
%
% the simplified constraint here accounts for causality and battery
% limitations
for jj = 1:H
    % simplified constraint
    combinedConstr(1+(jj-1)*N:N*jj) = energySpentPerSensorPerTimeStep(:,jj)...
                                     ...% all energy must go through battery: - sensorData.energyHarvestForOptimz(:,jj)... % current time step energy harvest
                                      - min(sensorData.Emax...
                                                ,...
                                            sensorData.energyRem... % the battery initial state
                                            + sensorData.energyHarvestForOptimz(:,jj)...
                                            + sum(EngyRcvdMinusSpent(:,1:jj-1),2)... % remainder of energy from prev time steps
                                            );

	sensorData.batteryStateOptmz(:,jj) = min(sensorData.Emax...
                                             ,...
                                             sensorData.energyRem... % the battery initial state
                                             + sum(EngyRcvdMinusSpent(:,1:jj),2)... % remainder of energy from prev time steps
                                            );
end

c = [totPowConstr;...
     combinedConstr;...
     -q]; % positivity of number of dynamic levels

ceq = [];

if sensorData.plotOptmzProgress
    recordProgress(totPowConstr,combinedConstr,q);
end

% =========================================================================
function recordProgress(totPowConstr,combinedConstr,q)
global progress

N = progress.numSensors;
H = progress.horizon;
numVarMod = length(q)/(N*H);


if isempty(progress.quantzLevels)
    progress.obj(1) =  abs(progress.objFunc(q));
    qvec = reshape(sum(reshape(q,numVarMod,N*H),1),N,H);
    progress.quantzLevels(:,1) = qvec(:);
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
    qvec = reshape(sum(reshape(q,numVarMod,N*H),1),N,H);
    progress.quantzLevels(:,end+1) = qvec(:);
    progress.constrValue(:,end+1) = [norm(totPowConstr);...
        min(abs(totPowConstr));...
        norm(combinedConstr);...
        min(abs(combinedConstr));...
        norm(-q);...
        ];
    progress.constrViol(:,end+1) = [1 * any(totPowConstr > 0);...
        2 * any(combinedConstr > 1e-8);...
        3 * any(-q > 1e-8)...
        ];
end
if progress.onlinePlotting && mod(size(progress.quantzLevels,2),200)==0
    ColorOrder = progress.ColorOrder;

    numIter = length(progress.obj);
    % draw the objective and norm of the variables
    figure(11)
    subplot(2,1,1)
    semilogy(1:numIter, abs(progress.obj),'-k+');hold on;
    semilogy(progress.quantzLevels','-');%,'Color',ColorOrder(1+mod(nn-1,21),:));
    hold off;
    subplot(2,1,2)
    plot(1:numIter, abs(progress.obj),'-k+');hold on;
    plot(progress.quantzLevels','-');%,'Color',ColorOrder(1+mod(nn-1,21),:));
    hold off;
    legend('MSE obj','Quantz Levels per node')
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
