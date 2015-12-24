function outputDataStruct = makeOutputDataStructEH(numSensors,dim,time,numConstraints,numStateVar)
% Creates a vector struct of sensor Kalman filter information.
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
% By Nick Roseveare, 03-2012

numElements = dim*2;
% outputDataStruct = struct('P',zeros(numElements,numElements),'Ptrace',[],'x',zeros(numElements,1),'t',[],'xEstHist',[],'optMSE',[],'actualMSE',[]);

% % outputDataStruct = repmat(outputDataStruct,numSensors,1);

tlength = length(time);
%outputDataStruct.optWeight      = zeros(numSensors*numStateVar,tlength);
outputDataStruct.x              = zeros(numElements,numSensors,tlength);
outputDataStruct.P              = inf*ones(2*dim,2*dim*numSensors,tlength);
outputDataStruct.Pdiag          = inf*ones(2*dim,numSensors,tlength);
outputDataStruct.Pmetric        = inf*ones(numSensors*numStateVar,numSensors,tlength);
outputDataStruct.xEstHist       = inf*ones(numElements,tlength);
outputDataStruct.xError         = inf*ones(1,tlength);
outputDataStruct.objectiveVal   = inf*ones(1,tlength);
% outputDataStruct.OptObjCompare  = inf*ones(1,tlength);
outputDataStruct.exitCode       = inf*ones(1,tlength);
outputDataStruct.energyUsed     = zeros(numSensors*numStateVar,tlength);
outputDataStruct.energyUsedPerSensor  = zeros(numSensors,tlength);
outputDataStruct.bits                 = zeros(numSensors*numStateVar,tlength);
outputDataStruct.bitsPerSensor        = zeros(numSensors,tlength);
outputDataStruct.batteryState   = zeros(numSensors,tlength);
% outputDataStruct.indivSensorMSE = inf*ones(numSensors,tlength);
outputDataStruct.t              = time;
% outputDataStruct.lagrMult       = inf*ones(numConstraints,tlength);
% outputDataStruct.activeSensors  = [];
% outputDataStruct.activeIndex    = [];

