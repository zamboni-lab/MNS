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
function outMat = mns_scanResult2table(mnsScanResults, model, sortBy, modelType)
if nargin < 3 || isempty(sortBy)
    sortBy = 'rankMax'; %options 'rankMax', 'rankSum', 'rankMaxRankSum'
end

if nargin < 4 || isempty(modelType)
    modelType = 'KEGGS';
end
%% calc sum and max breaking reaction and the product of both (rankMaxRankSum)
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
nlRange = [-1 mnsScanResults.nlRange];
maxBreakingTemp = maxBreaking;
maxBreaking = nlRange(mnsScanResults.maxBreakingReaction+1);
if isfield(mnsScanResults, 'maxSumBreakingReaction')
    maxSumBreaking = mnsScanResults.maxSumBreakingReaction;
else
    maxSumBreaking = sumBreaking.*maxBreaking;
end



%% get the tiedRanks for all reactions dependent on the sortBy value
isRankproduct = false;
switch sortBy
    case 'rankMax'
        breakingMetric = maxBreaking;
        rankBreaking = tiedrank(-breakingMetric);
    case 'rankSum'
        breakingMetric = sumBreaking;
        rankBreaking = tiedrank(-breakingMetric);
    case 'rankMaxRankSum'
        breakingMetric = maxSumBreaking;
        rankBreaking = tiedrank(-breakingMetric);
%     case 'rankDiff'
%         rankBreaking = tiedrank(-absDiffSubProdReaction);
%     case 'rankMaxRankDiff'
%         rankDiff = tiedrank(-absDiffSubProdReaction);
%         breakingMetric = maxBreaking;
%         rankBreaking = tiedrank(-breakingMetric);
%         rankBreaking = tiedrank(rankBreaking.*rankDiff);
%     case 'rankSumRankDiff'
%         rankDiff = tiedrank(-absDiffSubProdReaction);
%         breakingMetric = sumBreaking;
%         rankBreaking = tiedrank(-breakingMetric);
%         rankBreaking = tiedrank(rankBreaking.*rankDiff);
%     case 'rankMaxRankSumRankDiff'
%         rankDiff = tiedrank(-absDiffSubProdReaction);
%         breakingMetric = maxSumBreaking;
%         rankBreaking = tiedrank(-breakingMetric);
%         rankBreaking = tiedrank(rankBreaking.*rankDiff);
    case 'rankproductMax'
        rankBreaking = tiedrank(mnsScanResults.rankproductMaxBreaking);
        isRankproduct = true;
    case 'rankproductSum'
        rankBreaking = tiedrank(mnsScanResults.rankproductSumBreaking);
        isRankproduct = true;
    otherwise
        disp(['' sortBy '''is not a valid option. Results will be sorted by rankMax']);
        breakingMetric = maxBreaking;
        rankBreaking = tiedrank(-breakingMetric);
end
maxBreaking(maxBreaking == -1) = NaN;
%% generate the outtable
outHeader = [];
mns = mnsScanResults.mns{1};
outMat = [];

for i = 1:size(mns.results.reactionPairs,1)
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
%                     disp(idxEC(j))
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
        if isfield(model, 'merged') && model.merged == true
            idxRpMerged = find(model.mergedRp2mergedMetabolite(:,mns.results.reactionPairs(i,1)) ~= 0 ...
                | model.mergedRp2mergedMetabolite(:,mns.results.reactionPairs(i,2)) ~= 0);
            for r = 1:length(idxRpMerged)
                rpString = [rpString ', M:' model.rpId{idxRpMerged(r)}];
                idxEC = find(model.ECtoRP(:,idxRpMerged(r)) ~= 0);
                for j = 1:length(idxEC)
                    if ~isempty(idxEC)
                        %                     disp(idxEC(j))
                        idxGene = find(model.GeneToEC(:,idxEC(j)) ~= 0);
                    else
                        idxGene = [];
                    end

                    ECString = [ECString '; M: ' model.EC{idxEC(j)}];
                    for k = 1:length(idxGene)
                        geneSymbolString = [geneSymbolString '; M:' model.GeneSymbol{idxGene(k)}];
                    end
                end
                
            end
        end
        ECString = {ECString};
        rpString = {rpString};
        geneSymbolString = {geneSymbolString};
        switch modelType
            case 'KEGGS'
                if ~isRankproduct
                    if isfield(mnsScanResults, 'p')
                        switch mnsScanResults.p.pMode
                            case 'scan mode'
                                if sumBreaking(i) > 0
                                    pSum = mnsScanResults.p.pvalueSumBreaking(sumBreaking(i));
                                elseif sumBreaking(i) < 0
                                    pSum = mnsScanResults.p.pvalueSumBreaking(i);
                                else
                                    pSum = 1;
                                end
                                if maxBreakingTemp(i) > 0
                                    pMax = mnsScanResults.p.pvalueMaxBreaking(maxBreakingTemp(i));
                                elseif maxBreakingTemp(i) < 0
                                    pMax = mnsScanResults.p.pvalueMaxBreaking(i);
                                else
                                    pMax = 1;
                                end
                                if isfield(mnsScanResults.p, 'pvalueMaxSumBreaking')
                                    if maxSumBreaking(i) > 0
                                        pMaxSum = mnsScanResults.p.pvalueMaxSumBreaking(maxSumBreaking(i));
                                    else
                                        pMaxSum = 1;
                                    end
                                end
                                if isempty(outHeader)
                                    if isfield(mnsScanResults.p, 'pvalueMaxSumBreaking')
                                        outHeader = {'Rank' 'RP ID' 'EC' 'Gene Symbol' ...
                                            'Reaction' 'Max lambda1' 'p(max lambda1)' ...
                                            '#fractures' 'p(#fractures)' 'Max*Total Fractures' 'p(Max*Total Fractures)'};
                                    else
                                        outHeader = {'Rank' 'RP ID' 'EC' 'Gene Symbol' ...
                                            'Reaction' 'Max lambda1' 'p(max max lambda1)' ...
                                            '#fractures' 'p(#fractures)' 'Max*Total Fractures'};
                                    end
                                end
                                temp = [num2cell(rankBreaking(i)) rpString...
                                    ECString geneSymbolString ...
                                    {[model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}]}...
                                    num2cell(maxBreaking(i)) num2cell(pMax)...
                                    num2cell(sumBreaking(i)) num2cell(pSum) num2cell(maxSumBreaking(i))];
                                if isfield(mnsScanResults.p, 'pvalueMaxSumBreaking')
                                    temp = [temp num2cell(pMaxSum)];
                                end
                                outMat = [outMat; temp];
                            otherwise
                        end
                    else
                        if isempty(outHeader)
                            outHeader = {'Rank' 'RP ID' 'EC' 'Gene Symbol' ...
                                'Reaction' 'Max lambda1' 'Total Fractures' ...
                                'Max*Total Fractures'};
                        end
                        outMat = [outMat; num2cell(rankBreaking(i)) rpString...
                            ECString geneSymbolString ...
                            {[model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}]}...
                            num2cell(maxBreaking(i)) num2cell(sumBreaking(i)) num2cell(maxSumBreaking(i))];
                    end
                else
                    outHeader = {'Rank' 'RP ID' 'EC' 'Gene Symbol' ...
                        'Reaction' 'rankproduct(max lambda1)' 'p(rankproduct(max lambda1))' 'rankproduct(#fractures)' 'p(rankproduct(#fractures))'};
                    for j = 1:length(mnsScanResults.initMnsStruct)
                        outHeader{9+j} = ['rank(max lambda1) - para' num2str(j)];
                        outHeader{9+j+length(mnsScanResults.initMnsStruct)} = ['rank(#fractures) - para' num2str(j)];
                    end
                    outMat = [outMat; num2cell(rankBreaking(i)) rpString...
                            ECString geneSymbolString ...
                            {[model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}]}...
                            num2cell([mnsScanResults.rankproductMaxBreaking(i) mnsScanResults.pRankproductMaxBreaking(i) mnsScanResults.rankproductSumBreaking(i) ...
                            mnsScanResults.pRankproductSumBreaking(i) mnsScanResults.rankMaxBreaking(i,:) mnsScanResults.rankSumBreaking(i,:)] )];
                end
            otherwise
                disp('Not yet implemented');
                return;
        end
    end
end

[a,idxRankSort] = sort(cell2mat(outMat(:,1)));

outMat = [outHeader; outMat(idxRankSort,:)];

end