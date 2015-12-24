function state = filterMeasurement(state,z,Pz,t)
% This function filters a measurement update into the provided state.
%
% USAGE:
% state = filterMeasurement(state,z,Pz,t)
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

H = state.H;
P = state.P;
x = state.x;
statesize = length(x);

if isfield(state,'F')
    % This means that the time steps are uniform and a single F and Q are
    % adequate and do not need to be regenerated each time step
    F = state.F;
    Q = state.Q;
else
    tdiff = t - state.t;
    [F,Q] = makeFQ(2,tdiff,state.q,statesize/2 > 2);
end

% If the update time is sufficiently newer than the previous point, predict
% the current state
if tdiff > 1e-8
    % Update state to current time
    x = F*x;       % x[n+1|n]
    P = F*P*F' + Q;% P[n+1|n]
end

% Calculate kalman gain
Pzhat = H*P*H';
S = Pzhat + Pz;
% try
    L = chol(S,'lower'); % cholesky factorize because matrix Pos Semi def, cheaper "inverse" (i.e. back solve)
% catch
%     le = lasterror;
%     warning('MATLAB:filterMeasurement','diagonal loading necessary to stabilize covariance matrix: - last error was:\n %s',le.message)
%     ev = eig(P);
%     Pzhat = H*(P+eye(size(P))*(.01+abs(min(ev))))*H';
%     S = Pzhat + Pz;
%     L = chol(S+eye(size(S))*.01,'lower'); % cholesky factorize because matrix Pos Semi def, cheaper "inverse" (i.e. back solve)
% end
d = P*H'/L;
K = d/L';% kalman gain

% Measurement update
zhat = H*x;
state.x = x + K*(z-zhat);   % x[n|n] = x[n|n-1] + K[n]* z_innovation
state.P = P-K*H*P;          % P[n|n] = P[n|n-1] + K[n]* H * P
state.t = t;
if ~isfield(state,'recordHist') || state.recordHist == 1
    state.xHist(:,end+1) = state.x;
    state.PHist(:,:,end+1) = state.P;
    state.RHist(:,:,end+1) = Pz;
    state.zHist(:,end+1) = z;
    state.tHist(:,end+1) = t;
end

if isfield(state,'plotExtras')
    global covhist %#ok
    if ~isfield(covhist,'condP')
        covhist.condP = cond(P);
        covhist.trP = trace(P);
        covhist.condS = cond(S);
        covhist.innovNorm = norm(z-zhat);
    else
        rem = min(length(covhist.condP),148);%2.5*state.plotExtras.covhistLength * state.plotExtras.numActive);
        right = length(covhist.condP);
        left = max(1,right-rem);
        covhist.condP = [covhist.condP(left:right) cond(P)];
        covhist.trP = [covhist.trP(left:right) trace(P)];
        covhist.condS = [covhist.condS(left:right) cond(S)];
        covhist.innovNorm = [covhist.innovNorm(left:right) norm(z-zhat)];
    end
    if state.plotExtras.plotLog
        traceData = log10(covhist.trP);
        condData = log10(covhist.condP);
        innovNormData = log10(covhist.innovNorm);
        innovCovCond = log10(covhist.condS);
    else
        traceData = covhist.trP;
        condData = covhist.condP;
        innovNormData = covhist.innovNorm;
        innovCovCond = covhist.condS;
    end
    set(state.plotExtras.innovPlotHandle,'NextPlot','replaceChildren')
    set(state.plotExtras.innovPlotHandle,'GridLineStyle',':');
    plot(state.plotExtras.innovPlotHandle,innovNormData,'r','LineWidth',2)
    axis(state.plotExtras.innovPlotHandle,'auto')
    
    set(state.plotExtras.condPlotHandle,'NextPlot','replaceChildren')
    set(state.plotExtras.innovPlotHandle,'GridLineStyle',':');
    plot(state.plotExtras.condPlotHandle,condData,'k','LineWidth',1.5);hold on;
    set(state.plotExtras.condPlotHandle,'NextPlot','add')
    plot(state.plotExtras.condPlotHandle,traceData,'b--','LineWidth',1.5);hold on;
    
    plot(state.plotExtras.condPlotHandle,innovCovCond,'r','LineWidth',1.5);hold off;
    axis(state.plotExtras.condPlotHandle,'auto')
end



