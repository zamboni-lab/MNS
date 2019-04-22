%% Copyright and License
% Metabolic Network Segmentation Toolbox – A probabilistic graphical 
% modelling tool to identify sites and sequential order of metabolic regulations
% Copyright (C) 2016, Andreas Kühne & Nicola Zamboni
% 
% This file is part of the Metabolic network segmentation toolbox
% 
% Metabolic network segmentation toolbox is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation, either
% version 3 of the License or any later version.
% 
% Metabolic network segmentation toolbox is distributed in the hope 
% that it will be useful, but WITHOUT ANY WARRANTY; without even the
% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy (gpl.txt) of the GNU General Public License
% along with the metabolic network segmentation toolbox. If not,
% see http://www.gnu.org/licenses/.
%
function results = mns_loadResult2state(folder)
fRes = 'out.txt';
fSolution = 'outSolution.txt';
fVarMarg = 'outMarginals.txt';
fFacMarg = 'outFactorMarginals.txt';

marg = false;

%open Solution File
fid = fopen(fullfile(folder, fRes),'r');
res = fscanf(fid, '%d\t%d');
results.varId = res(1:2:size(res,1)-1);
results.inferedLabels = res(2:2:size(res,1));
fclose(fid);
clear res

% get solution value and bound
fid = fopen(fullfile(folder, fSolution),'r');
line = fgetl(fid);
while(line ~= -1)
    if regexp(line,'Solution')
       solution = regexprep(line, 'Solution', '');
       solution = regexprep(solution, '\s*\t*', '');
       results.solution = str2num(solution);
    elseif regexp(line,'Bound')
        solution = regexprep(line, 'Bound', '');
       solution = regexprep(solution, '\s*\t*', '');
       if regexp(solution,'INF')
           results.bound = inf;
       else
           results.bound = str2num(solution);
       end
       
    end
    line = fgetl(fid);
end
fclose(fid);

%get variable marginals
fid = fopen([folder fVarMarg],'r');
if fid ~= -1
    res = fscanf(fid, '%d\t%d\t%f');
    varid = res(1:3:size(res,1)-2);
    results.labelId = unique(res(2:3:size(res,1)-1));
    mLab = max(results.labelId);
    marginals = res(3:3:size(res,1));
    results.varMarginals = reshape(marginals,mLab+1,max(varid)+1)';
    
    clear res
    fclose(fid);
end


%get factor marginals
fid = fopen([folder fFacMarg],'r');
if fid ~= -1
    res = fscanf(fid, '%d\t%d\t%f');
    facid = res(1:3:size(res,1)-2);
    configid = res(2:3:size(res,1)-1);
    facMarginals = res(3:3:size(res,1));
    results.facId = unique(facid);
    results.configId = unique(configid);
    results.facMarginalsMat = nan(length(results.facId),length(results.configId));
    for i = 1:length(results.facId)
        results.facMarginals(i).facId = results.facId(i);
        results.facMarginals(i).configid = configid(facid == results.facId(i));
        results.facMarginals(i).facPotentials = facMarginals(facid == results.facId(i));
        results.facMarginalsMat(i,1:length(find(facid == results.facId(i)))) = results.facMarginals(i).facPotentials;
    end
    fclose(fid);
end



end