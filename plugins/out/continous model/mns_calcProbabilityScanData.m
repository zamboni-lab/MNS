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
function [probArr, score, idx, idxScore, idxVec] = mns_calcProbabilityScanData(mnsResults, model,ws,wn,plotDistributions, plotMnsResults)
if nargin < 3
    ws = 0.1;
    wn = 0.25;
end
if nargin < 5
    plotDistributions = true;
end
if nargin < 6 
    plotMnsResults = true;
end
    
%%
uNL1 = unique(mnsResults.nL1);
uTL3 = unique(mnsResults.tL3);
prob = zeros(length(mnsResults.mns),1);
probArr = zeros(length(uNL1),length(uTL3));

% tic
% for i = 1:length(mnsResults.mns)
%     for j = 1:length(mnsResults.mns{i}.results.inferedLabels)
%         idObsClique = find(mnsResults.mns{i}.observations.modelVarId == j-1);
%         prob(i) = prob(i)+ mnsResults.mns{i}.observations.obsCliquePotential(idObsClique,mnsResults.mns{i}.results.inferedLabels(j)+1);
%     end
%     %
%     probArr(uNL1 == mnsResults.nL1(i), uTL3 == mnsResults.tL3(i)) = prob(i);
% end
% toc
% tic
idObsClique = zeros(length(mnsResults.mns{1}.results.inferedLabels),1);
for j = 1:length(mnsResults.mns{1}.results.inferedLabels)
    idObsClique(j,1) = find(mnsResults.mns{1}.observations.modelVarId == j-1);
end
for i = 1:size(mnsResults.clustLabel,3)
%     linearInd = sub2ind(size(mnsResults.mns{i}.observations.obsCliquePotential), idObsClique, mnsResults.mns{i}.results.inferedLabels+1);
    inferedLabels = reshape(mnsResults.clustLabel(:,:,i),length(mnsResults.mns{1}.results.inferedLabels),1);
    linearInd = sub2ind(size(mnsResults.mns{1}.observations.obsCliquePotential), idObsClique, inferedLabels+1);
    probArr(uNL1 == mnsResults.nL1(i), uTL3 == mnsResults.tL3(i)) = sum(mnsResults.mns{1}.observations.obsCliquePotential(linearInd));
%     probArr(uNL1 == mnsResults.nL1(i), uTL3 == mnsResults.tL3(i)) = sum(mnsResults.mns{i}.observations.obsCliquePotential(linearInd));
end
% toc
%% calculate score
% normalize probArr
probArr = probArr-(min(min(probArr)));
% old 
wt = 1/max(max(mnsResults.amountTemporalBreaks))*ws;
wf = 1/max(max(mnsResults.amountBreakingReaction))*wn;
wo = 1/max(max(probArr));

% new
maxTempBreaks = size(mnsResults.temporalBreaks,1)*size(mnsResults.temporalBreaks,2);
wt = 1/maxTempBreaks*ws;
maxBreaks = size(mnsResults.breakingReaction,1)*size(mnsResults.breakingReaction,2);
wf = 1/maxBreaks*wn;

wo = 1/max(max(probArr));
score = probArr*wo - mnsResults.amountTemporalBreaks*wt - mnsResults.amountBreakingReaction*wf;
%change: we don't want to have fittings that don't explain the data at
%all?!
score(probArr == 0 | mnsResults.amountTemporalBreaks == 0 | mnsResults.amountBreakingReaction == 0) = NaN;
if plotDistributions
    mns_plotTemporalScreeningResults(mnsResults,probArr, score,false)
end
%%
idxVec = find(score == max(max(score)));
% idxVec = idx;
% idx = idx(1);
% maxTL3 = ceil(idx/size(score,1));
% maxNl1 = mod(idx, size(score,1));
% if maxNl1 == 0
%     maxNl1 = size(score,1);
% end

for i = 1:length(idxVec)
    [maxNl1, maxTL3] = ind2sub(size(score),idxVec(i));
    idxScore = [ maxNl1, maxTL3];
    idxVec(i) = find(mnsResults.nL1 == uNL1(maxNl1) & mnsResults.tL3 == uTL3(maxTL3));
end
idx = idxVec(1);
if plotMnsResults
    mns_plotTemporalDataResults(mnsResults,model,idx, sprintf('ws=%0.3f, wn = %0.3f',ws,wn))
end

end
