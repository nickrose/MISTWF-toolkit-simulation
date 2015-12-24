function sample = fixedSampleDatabase(dbid,rvGenFunc,reuse,idx,idx2,data)

% fixedSampleDatabase.m - this function stores a database of (noise or
% other) samples for the purpose of keeping noise samples exactly the same
% across different formulations of the problem so that solutions are
% comparable.
% NOTE: This function relies on the fact that a (noise) sample for any
% particular distribution should always have instance value |s| < inf.
% This EXCLUDES the case where the 'data' option is being used.
%
% USAGE:
%     sample = fixedSampleDatabase(dbid,rvGenFunc,reuse,idx,idx2,data)
%
% INPUTS:
%  dbid         input string identifying which database of noise samples is
%               to be accessed/modified
%  rvGenFunc    the function handle used to generate the random variable
%               (e.g. @randn, @rand, etc.), can also be 'data' which means
%               the data input should be utilized and (when reuse is 1)
%               will be inserting data instead of drawing random values
%  reuse        specifies whether the database entry should be retrieved or
%               regenerated
%  idx          the column index into the specified database (this is
%               usually the only index specified, but there is any option
%               to keep a two dimensional database), it is allowable to
%               specify an array of indices into the columns
%  idx2         (optional) specify the row index for the 2-D array option,
%               we must assume that this index is always a singleton for
%               proper operation, generation, and lookup
%  data         (optional) the data necessary to insert in the case when
%               rvGenFunc is 'data', WARNING: length of data should match
%               the length of idx
%
% OUTPUT:
%
%  sample       the retrieved/generated sample value
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
% By: Nick Roseveare, Feb 2011

persistent sampleDB
if nargin == 1
    if strcmp(dbid,'cleardb')
        sampleDB = [];
    end
    if nargout > 0
        sample = 'dbcleared';
    end
else
    if nargin < 5 || (strcmp(rvGenFunc,'data') && ((reuse ~= 1 && nargin < 7) || (reuse == 1 && nargin < 6)))
        % For the 1-D array case, just let this always be 1
        idx2 = 1;
    end
    
    % If the database is empty, create the first entry
    if isempty(sampleDB)
        sampleDB = struct(dbid,[inf]);
    elseif ~isfield(sampleDB,dbid)
        sampleDB.(dbid) = inf;
    end
    
    numInd = length(idx);
    [currDBdepth, currDBlen] = size(sampleDB.(dbid));
    iterDBlen = currDBlen;
    if idx2 > currDBdepth
        sampleDB.(dbid)(currDBdepth+1:idx2,:) = inf*ones(idx2-currDBdepth,currDBlen);
    end
    
    % Check to see if the given index is beyond the size of the database OR if
    % it has not been specified OR if the values are not being recycled
    if ~strcmp(rvGenFunc,'data')
        if currDBlen < max(idx) || ... % <- this condition is what keeps the if statement, otherwise we would just use the update indices vector
                any(sampleDB.(dbid)(idx2,idx(idx <= currDBlen)) == repmat(inf,1,sum(idx <= currDBlen))) || ...
                ~reuse
            updatedIndices = idx(idx > currDBlen ...
                | [(sampleDB.(dbid)(idx2,idx(idx <= currDBlen)) == repmat(inf,1,sum(idx <= currDBlen))) repmat(1,1,sum(idx > currDBlen))] ...
                | repmat(reuse,1,numInd) == 0);
            for ii = updatedIndices
                % Fill in skipped locations with inf (i.e. representing unallocated)
                sampleDB.(dbid)(:,iterDBlen+1:ii) = inf;
                iterDBlen = size(sampleDB.(dbid),2);
                
                % Generate the sample
                sampleDB.(dbid)(idx2,ii) = rvGenFunc(1,1);
            end
        end
    else
        % Changed data entering so that even if 'reuse' is on, will fill in the
        % database if it is empty
        if currDBlen < max(idx) || ... % <- this condition is what keeps the if statement, otherwise we would just use the update indices vector
                any(sampleDB.(dbid)(idx2,idx(idx <= currDBlen)) == repmat(inf,1,sum(idx <= currDBlen))) || ...
                ~reuse % store values
            
            updatedIndices = idx(idx > currDBlen ...
                | [(sampleDB.(dbid)(idx2,idx(idx <= currDBlen)) == repmat(inf,1,sum(idx <= currDBlen))) repmat(1,1,sum(idx > currDBlen))] ...
                | repmat(reuse,1,numInd) == 0);
            for ii = updatedIndices
                % Fill in skipped locations with inf (i.e. representing unallocated)
                sampleDB.(dbid)(:,iterDBlen+1:ii) = inf;
                iterDBlen = size(sampleDB.(dbid),2);
                
                % Fill in with provided data - should not be adjusting more
                sampleDB.(dbid)(idx2,ii) = data(ii == idx);
            end
        else
            % test to make sure indices are valid
            if numel(sampleDB.(dbid)) < max(idx)*idx2
                error('The indices provide are for data which is unassigned in the sample database')
            end
        end
    end
    
    % Return the sample(s)
    if nargout == 1
        sample = sampleDB.(dbid)(idx2,idx);
    end
end