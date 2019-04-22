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
function mns_scanResultHits2cytoscape(model, mnsScanResults, nameBase, sortBy, idx, nNeighborMetabolites)
if nargin < 3 || isempty(nameBase)
    nameBase = 'temp_mns';
end
if nargin < 4 || isempty(sortBy)
    sortBy = 'rankMax';
end
if nargin < 5 || isempty(idx)
    idx = 1;
end
if nargin < 6 || isempty(nNeighborMetabolites)
    nNeighborMetabolites = 10;
end
%%
sPaths = graphallshortestpaths(sparse(model.iaMat), 'Directed' ,0);
%% get sorted results
outTable = mns_scanResult2table(mnsScanResults, model,sortBy);
rpIds = outTable(idx+1,2);
geneSymbols = outTable(idx+1,4);
ranks = outTable(idx+1,1);
metIdx = [];
edgeLabelArr = cell(length(idx), 1);

for i = 1:length(rpIds)
    rpModelIdx(i) = find(strcmp(model.rpId, rpIds{i}) == 1);
    nMetIdx = find(model.mat(rpModelIdx(i),:) ~= 0);
    nMetIdxArr(i,:) = nMetIdx;
    edgeLabelArr{i} = [num2str(ranks{i}) ': ' geneSymbols{i}];
    sPathsTemp = min(sPaths(nMetIdx,:));
    temp = sort(sPathsTemp);
    distCo = temp(nNeighborMetabolites);
    if isinf(distCo)
        distCo = max(temp(isinf(temp) == 0));
    end
    metIdx = [metIdx find(sPathsTemp <= distCo)];
end
metIdx = unique(metIdx);
idxRp2plot = find(sum(model.mat(:,metIdx),2)>1);
%% Results 2 cytoscape
modelType = 'KEGGS';

edgeFile = [nameBase ' - mns_hits_edges.xls'];
nodeFile = [nameBase ' - mns_hits_nodes.xls'];

nMet = length(model.metaboliteId);
mns = mnsScanResults.mns{1};
metOutMat(:,1:4) = [{'MetId'} {'metName'} {'Detected'} {'Data'}; ...
    model.metaboliteId(metIdx) model.metaboliteName(metIdx) num2cell(mns.observations.detected(metIdx))...
    num2cell(mns.observations.data(metIdx))];
pairOutMat = {'Substrate' 'rpId' 'Product' , 'EdgeLabel', 'EC', 'GeneSymbol'...
    'Rank(TotalFractures)' 'TotalFractures' 'Rank(MaxFractures)' ...
    'MaxFractures' 'Rank(Max*TotalFractures)' 'Rank(Delta(Sub,Prod))' 'Delta(Sub,Prod)'};


mnsRpTemp = mnsScanResults.mns{1}.results.reactionPairs;


sumBreaking = sum(mnsScanResults.breakingReaction,2);
maxBreaking = zeros(size(sumBreaking));
for i = 1:length(maxBreaking)
    idx = find(mnsScanResults.breakingReaction(i,:) == 1);
    if ~isempty(idx)
        maxBreaking(i) = max(idx(end));
    end
end
rankSumBreaking = tiedrank(-sumBreaking);
rankMaxBreaking = tiedrank(-maxBreaking);

if isfield(mnsScanResults, 'maxSumBreakingReaction')
    maxSumBreaking = mnsScanResults.maxSumBreakingReaction;
else
    maxSumBreaking = sumBreaking.*maxBreaking;
end

if isfield(mnsScanResults, 'absDiffSubProdReaction')
    absDiffSubProdReaction = mnsScanResults.absDiffSubProdReaction;
else
    absDiffSubProdReaction = zeros(length(sumBreaking),1);
end

rankMaxSumBreaking = tiedrank(-maxSumBreaking);
rankAbsDiffSubProdReaction = tiedrank(-absDiffSubProdReaction);
edgeLabel = cell(length(idxRp2plot),1);
%%
c = 1;
for i = 1:size(idxRp2plot)
    if isfield(model, 'mat')
%         idxRp = find(model.mat(:,mns.results.reactionPairs(i,1)) ~= 0 & model.mat(:,mns.results.reactionPairs(i,2)) ~= 0);
        idxRp = idxRp2plot(i);
        idxMets = sort(find(model.mat(idxRp,:) ~= 0));
        
        idxMns = find((mns.results.reactionPairs(:,1) == idxMets(1) &  mns.results.reactionPairs(:,2) == idxMets(2)) ...
            | (mns.results.reactionPairs(:,1) == idxMets(2) &  mns.results.reactionPairs(:,2) == idxMets(1)));
        
        ECString = '';
        geneSymbolString = '';
        idxRankedRp = find(rpModelIdx == idxRp);
        if isempty(idxRankedRp)
            edgeLabel{i,1} = '';
        else
            edgeLabel{i,1} = edgeLabelArr{idxRankedRp};
        end
        
        for r = 1:length(idxRp)
            if r == 1
                rpString = model.rpId{idxRp(r)};
            else
                rpString = [rpString ', ' model.rpId{idxRp(r)}];
            end
            idxEC = find(model.ECtoRP(:,idxRp(r)) ~= 0);
            for j = 1:length(idxEC)
                if ~isempty(idxEC)
                    idxGene = find(model.GeneToEC(:,idxEC(j)) ~= 0);
                else
                    idxGene = [];
                end
                if j == 1 && r == 1
                    ECString = model.EC{idxEC(j)};
                else
                    ECString = [ECString '; ' model.EC{idxEC(j)}];
                end
                for k = 1:length(idxGene)
                    if k == 1 && j == 1 && r == 1
                        geneSymbolString = model.GeneSymbol{idxGene(k)};
                    else
                        geneSymbolString = [geneSymbolString '; ' model.GeneSymbol{idxGene(k)}];
                    end
                end
            end
            
        end
        ECString = {ECString};
        rpString = {rpString};
        geneSymbolString = {geneSymbolString};
        
    end
%     
%     if i == 107
%         disp(i)
%     end
    if ~isempty(idxMns)
        c = c+1;
        pairOutMat(c,:) = [model.metaboliteId(idxMets(1)) ...
            rpString model.metaboliteId(idxMets(2)) edgeLabel(i) ECString geneSymbolString ...
            num2cell(rankSumBreaking(idxMns)) num2cell(sumBreaking(idxMns)) ...
            num2cell(rankMaxBreaking(idxMns)) num2cell(maxBreaking(idxMns)) ...
            num2cell(rankMaxSumBreaking(idxMns)) num2cell(rankAbsDiffSubProdReaction(idxMns))...
            num2cell(absDiffSubProdReaction(idxMns))];
    end
end

% for i = 1:length(mnsScanResults.mns)
%     pairOutMat(:,end+1) = [{['nl_' num2str(i)]}; num2cell(mnsScanResults.mns{i}.results.breakingReactions(idxRp2plot))];
%     metOutMat(:,end+1) = [{['nl_' num2str(i)]}; num2cell(mnsScanResults.mns{i}.results.inferedLabels(metIdx))];
% end

xlswrite(edgeFile, pairOutMat);
xlswrite(nodeFile, metOutMat);


end