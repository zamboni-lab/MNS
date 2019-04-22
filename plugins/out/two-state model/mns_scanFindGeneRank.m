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
function [rankSum, rankMax, rankUSum, rankUMax, maxRankUsum, maxRankUMax, gene2geneDist, gene2rp] = mns_scanFindGeneRank(mnsScanResults, model,geneSymbol, k)
if isfield(mnsScanResults, 'combinedParameterizations') && mnsScanResults.combinedParameterizations
    isRankProduct = true;
else
    isRankProduct = false;
end

if ~isfield(model, 'gene2geneDist') || ~isfield(model, 'gene2rp')
    disp('calc matrix');
    [model.gene2geneDist,model.gene2rp] = mns_gene2geneDistMat(model);
end
gene2geneDist = model.gene2geneDist; 
gene2rp = model.gene2rp;

%calc sum and max breaking reaction
if isfield(mnsScanResults, 'sumBreakingReaction')
    sumBreaking = mnsScanResults.sumBreakingReaction;
else
    sumBreaking = sum(mnsScanResults.breakingReaction,2);
end

if isfield(mnsScanResults, 'maxBreakingReaction')
    maxBreaking = mnsScanResults.maxBreakingReaction;
else
    maxBreaking = zeros(size(sumBreaking));
    for i = 1:length(maxBreaking)
        idx = find(mnsScanResults.breakingReaction(i,:) == 1);
        if ~isempty(idx)
            maxBreaking(i) = max(idx(end));
        end
    end
end
if ~isRankProduct
    rankSumBreaking = tiedrank(-sumBreaking);
    rankMaxBreaking = tiedrank(-maxBreaking);
    uSumBreaking = unique(sumBreaking);
    rankUSumBreaking = tiedrank(-uSumBreaking);
    uMaxBreaking = unique(maxBreaking);
    rankUMaxBreaking = tiedrank(-uMaxBreaking);
else
    rankSumBreaking = tiedrank(mnsScanResults.rankproductSumBreaking);
    rankMaxBreaking = tiedrank(mnsScanResults.rankproductMaxBreaking);
    sumBreaking = mnsScanResults.rankproductSumBreaking;
    uSumBreaking = unique(sumBreaking);
    rankUSumBreaking = tiedrank(uSumBreaking);
    maxBreaking = mnsScanResults.rankproductMaxBreaking;
    uMaxBreaking = unique(maxBreaking);
    rankUMaxBreaking = tiedrank(uMaxBreaking);
end
maxRankUsum = rankUSumBreaking(1);
maxRankUMax = rankUMaxBreaking(1);
rankSumBreakingRp = inf(length(model.rpId),1);
rankMaxBreakingRp = inf(length(model.rpId),1);
rankUSumBreakingRp = inf(length(model.rpId),1);
rankUMaxBreakingRp = inf(length(model.rpId),1);

mns = mnsScanResults.mns{1};
for i = 1:size(mns.results.reactionPairs,1)
    if isfield(model, 'mat')
        idxRp = find(model.mat(:,mns.results.reactionPairs(i,1)) ~= 0 & model.mat(:,mns.results.reactionPairs(i,2)) ~= 0);
        rankSumBreakingRp(idxRp) = rankSumBreaking(i);
        rankMaxBreakingRp(idxRp) = rankMaxBreaking(i);
        rankUSumBreakingRp(idxRp) = rankUSumBreaking(uSumBreaking == sumBreaking(i));
        rankUMaxBreakingRp(idxRp) = rankUMaxBreaking(uMaxBreaking == maxBreaking(i));
    end
end

gene2geneTemp = model.gene2geneDist((strcmp(model.GeneSymbol,geneSymbol) == 1),:);
rankSum = inf(1,k+1);
rankMax = inf(1,k+1);
rankUSum = inf(1,k+1);
rankUMax = inf(1,k+1);
for i = 0:k
    idxGene = find(gene2geneTemp <= i);
    for j = 1:length(idxGene)
        rankSum(i+1) = min([rankSum(i+1); rankSumBreakingRp(model.gene2rp(idxGene(j),:) == 1)]);
        rankMax(i+1) = min([rankMax(i+1); rankMaxBreakingRp(model.gene2rp(idxGene(j),:) == 1)]);
        rankUSum(i+1) = min([rankUSum(i+1); rankUSumBreakingRp(model.gene2rp(idxGene(j),:) == 1)]);
        rankUMax(i+1) = min([rankUMax(i+1); rankUMaxBreakingRp(model.gene2rp(idxGene(j),:) == 1)]);
    end
end


end