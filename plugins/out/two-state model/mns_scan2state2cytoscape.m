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
function mns_scan2state2cytoscape(model, mnsResults, nameBase)

modelType = 'KEGGS';

edgeFile = [nameBase ' - mns_edges.xls'];
nodeFile = [nameBase ' - mns_nodes.xls'];

nMet = length(model.metaboliteId);
mns = mnsResults.mns{1};
metOutMat(:,1:4) = [{'MetId'} {'metName'} {'Detected'} {'Data'}; ...
    model.metaboliteId model.metaboliteName num2cell(mns.observations.detected)...
    num2cell(mns.observations.data)];
pairOutMat = {'Substrate' 'rpId' 'Product' ,'EC', 'GeneSymbol'...
    'Rank(TotalFractures)' 'TotalFractures' 'Rank(MaxFractures)' ...
    'MaxFractures' 'Rank(Max*TotalFractures)' 'Rank(Delta(Sub,Prod))' 'Delta(Sub,Prod)'};


mnsRpTemp = mnsResults.mns{1}.results.reactionPairs;


sumBreaking = sum(mnsResults.breakingReaction,2);
maxBreaking = zeros(size(sumBreaking));
for i = 1:length(maxBreaking)
    idx = find(mnsResults.breakingReaction(i,:) == 1);
    if ~isempty(idx)
        maxBreaking(i) = max(idx(end));
    end
end
rankSumBreaking = tiedrank(-sumBreaking);
rankMaxBreaking = tiedrank(-maxBreaking);

if isfield(mnsResults, 'maxSumBreakingReaction')
    maxSumBreaking = mnsResults.maxSumBreakingReaction;
else
    maxSumBreaking = sumBreaking.*maxBreaking;
end

if isfield(mnsResults, 'absDiffSubProdReaction')
    absDiffSubProdReaction = mnsResults.absDiffSubProdReaction;
else
    absDiffSubProdReaction = zeros(length(sumBreaking),1);
end

rankMaxSumBreaking = tiedrank(-maxSumBreaking);
rankAbsDiffSubProdReaction = tiedrank(-absDiffSubProdReaction);

for i = 1:size(mnsResults.mns{1}.results.reactionPairs,1)
    if isfield(model, 'mat')
        idxRp = find(model.mat(:,mns.results.reactionPairs(i,1)) ~= 0 & model.mat(:,mns.results.reactionPairs(i,2)) ~= 0);
        
        ECString = '';
        geneSymbolString = '';
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
%     disp(i)
    pairOutMat(i+1,:) = [model.metaboliteId(mnsRpTemp(i,1)) ...
        rpString model.metaboliteId(mnsRpTemp(i,2)) ECString geneSymbolString ...
        num2cell(rankSumBreaking(i)) num2cell(sumBreaking(i)) ...
        num2cell(rankMaxBreaking(i)) num2cell(maxBreaking(i)) ...
        num2cell(rankMaxSumBreaking(i)) num2cell(rankAbsDiffSubProdReaction(i))...
        num2cell(absDiffSubProdReaction(i))]; 
end

for i = 1:length(mnsResults.mns)
    pairOutMat(:,end+1) = [{['nl_' num2str(i)]}; num2cell(mnsResults.mns{i}.results.breakingReactions)];
    metOutMat(:,end+1) = [{['nl_' num2str(i)]}; num2cell(mnsResults.mns{i}.results.inferedLabels)];
end

xlswrite(edgeFile, pairOutMat);
xlswrite(nodeFile, metOutMat);

end