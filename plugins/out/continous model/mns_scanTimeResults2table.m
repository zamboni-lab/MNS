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
function [seqFracList, seqFracListOnlyFracPos, nFracList, nFracListSplit] ...
    = mns_scanTimeResults2table(mnsResults,model, dWs,dWn,limit, mode)
if nargin < 3
    dWs = 0.01;
end
if nargin < 4
    dWn = 0.01;
end
if nargin < 5 || isempty(limit)
    limit = inf;
end
if nargin < 6
    mode = 'scan';
end

%% start scan
fractures = true;
wt = 0;
wn = 0;
c = 0;
nFrames = size(mnsResults.clustLabel,2);
seqFracList = []; %contains results with different frac direction and label span
seqFracListOnlyFracPos = []; %contains results w/o frac direction and label span
nFracList = []; %total number of fractures
nFracListSplit = []; %split number of fractures
while fractures
    wt = wt+dWs;
    wn = wn+dWn;
    
    disp(wt)
    c = c + 1;
    [prob, score, idx, idxScore] = mns_calcProbabilityScanData(mnsResults, model,wn,wt,false, false);
    
    if mnsResults.amountBreakingReaction(idxScore(1),idxScore(2)) == 0 ...
            && mnsResults.amountTemporalBreaks(idxScore(1),idxScore(2)) == 0
        fractures = false;
        break;
    end
    if wt >= limit || wn >= limit
        fractures = false;
        break;
    end
    
    % initialize frequencies
    freqLengthLabelSpan(c,:) = zeros(1,nFrames);
    freqTempFractures(c,:) = zeros(1,nFrames);
    freqNFractures(c,:) = zeros(1,nFrames+1);
    freqLengthNFrac(c,:) = zeros(1,nFrames+1);
    % get data for max score
    labels = mnsResults.clustLabel(:,:,idx);
    breakingReaction = mnsResults.breakingReaction(:,:,idx);
    temporalBreaks = mnsResults.temporalBreaks(:,:,idx);
    
    %% go through temporal breaks
    for i = 1:size(temporalBreaks,1) % go through all metabolites
        idxTempFrac = find(temporalBreaks(i,:) == 1);
        if ~isempty(idxTempFrac)
            freqTempFractures(c,length(idxTempFrac)+1) = freqTempFractures(c,length(idxTempFrac)+1)+1;
            for j = 1:length(idxTempFrac) % go through all breaks
                %%
                % calc difference between labels to indicate in which
                % direction the metabolite changes delta label > 0 increase
                % < 0 decrease
                deltaLabel = labels(i, idxTempFrac(j)+1)-labels(i, idxTempFrac(j));
                
                % length of label describes how the number of frames to the
                % previous and next temporal fracture
                
                if j == 1
                    lengthLabel = 0;
                else
                    lengthLabel = idxTempFrac(j-1);
                end
                if j == length(idxTempFrac)
                    lengthLabel = nFrames-lengthLabel;
                else
                    lengthLabel = idxTempFrac(j+1)-lengthLabel;
                end
                freqLengthLabelSpan(c,lengthLabel) = freqLengthLabelSpan(c,lengthLabel) + 1;
                tempFrac = [i idxTempFrac(j)+0.5 deltaLabel lengthLabel c wn wt];
                
                %% search and update list
                if c == 1
                    seqFracList = [seqFracList; tempFrac];
                    seqFracListOnlyFracPos = [seqFracListOnlyFracPos; tempFrac];
                else
                    idxRow = find(ismember(seqFracList(:,1:4),tempFrac(1:4),'rows'),1);
                    if ~isempty(idxRow)
                        seqFracList(idxRow,:) = tempFrac;
                    else
                        seqFracList = [seqFracList; tempFrac];
                    end
                    idxRow = find(ismember(seqFracListOnlyFracPos(:,1:2),tempFrac(1:2),'rows'),1);
                    if ~isempty(idxRow)
                        seqFracListOnlyFracPos(idxRow,:) = tempFrac;
                    else
                        seqFracListOnlyFracPos = [seqFracListOnlyFracPos; tempFrac];
                    end
                end
                
            end
        else
            freqLengthLabelSpan(c,1) = freqLengthLabelSpan(c,1) + 1;
            freqTempFractures(c,1) = freqTempFractures(c,1)+1;
        end
    end
    %calc probabilits to find reactions with >= labelspan or >= number
    %of fractures by chance
    for j = 1:nFrames
        pLengthLabelSpan(c,j) = sum(freqLengthLabelSpan(c,j:end))/sum(freqLengthLabelSpan(c,:));
        pFractures(c,j) = sum(freqTempFractures(c,j:end))/sum(freqTempFractures(c,:));
    end
    %% go through all neighborhood breaking reaction
    for i = 1:size(breakingReaction,1) % go through all reactions
        idxNFrac = find(breakingReaction(i,:) == 1);
        if ~isempty(idxNFrac)
%             disp(length(idxNFrac)+1)
            freqNFractures(c,length(idxNFrac)+1) = freqNFractures(c,length(idxNFrac)+1)+1;
            fracGr = zeros(length(idxNFrac),1);
            for j = 1:length(idxNFrac)
                if j == 1
                    fracGr(j) = 1;
                else
                    if idxNFrac(j)-1 == idxNFrac(j-1)
                        fracGr(j) = max(fracGr);
                    else
                        fracGr(j) = max(fracGr)+1;
                    end
                end
            end
            uFracGr = unique(fracGr);
%             nFracReac = [i length(idxNFrac) length(uFracGr) {idxNFrac} c wn wt];
            nFracReac = [i length(idxNFrac) length(uFracGr) c wn wt];
            %% search and update list
            if c == 1
                nFracList = [nFracList; nFracReac];
            else
                idxRow = find(ismember(nFracList(:,1:3),nFracReac(1:3),'rows'),1);
                if ~isempty(idxRow)
                    nFracList(idxRow,:) = nFracReac;
                else
                    nFracList = [nFracList; nFracReac];
                end
            end
            for j = 1:length(uFracGr)
                idxTemp = idxNFrac(uFracGr(j) == fracGr);
                lengthFrac = length(idxTemp);
                freqLengthNFrac(c,lengthFrac+1) = freqLengthNFrac(c,lengthFrac+1)+1;
                nFrac = [i min(idxTemp) max(idxTemp) lengthFrac c wn wt];
                %% search and update list
                if c == 1
                    nFracListSplit = [nFracListSplit; nFrac];
                else
                    idxRow = find(ismember(nFracListSplit(:,1:4),nFrac(1:4),'rows'),1);
                    if ~isempty(idxRow)
                        nFracListSplit(idxRow,:) = nFrac;
                    else
                        nFracListSplit = [nFracListSplit; nFrac];
                    end
                end
            end
        else
            freqNFractures(c,1) = freqNFractures(c,1)+1;
            freqLengthNFrac(c,1) = freqLengthNFrac(c,1)+1;
        end
    end
    %calc probabilits to find reactions with >= labelspan or >= number
    %of fractures by chance
    for j = 1:nFrames
        pLengthLabelSpan(c,j) = sum(freqLengthLabelSpan(c,j:end))/sum(freqLengthLabelSpan(c,:));
        pFractures(c,j) = sum(freqTempFractures(c,j:end))/sum(freqTempFractures(c,:));
    end
    switch mode
        case 'exact'
            break;
    end
end

%% add metabolite names etc
[~,idxTemp] = sort(seqFracList(:,5),'descend');
tempHeader = {'Rank' 'Metabolite Name' 'Metabolite Id' 'Metabolite Model ID' 'Fracture Position'...
    'module change' , 'Module Label Span', 'step', 'wn', 'wt'};
    
seqFracList = [tempHeader; num2cell(tiedrank(-seqFracList(idxTemp,5))), model.metaboliteName(seqFracList(idxTemp,1)),...
    model.metaboliteId(seqFracList(idxTemp,1)), num2cell(seqFracList(idxTemp,:))];
%%
[~,idxTemp] = sort(seqFracListOnlyFracPos(:,5),'descend');
tempHeader = {'Rank' 'Metabolite Name' 'Metabolite Id' 'Metabolite Model ID' 'Fracture Position'...
    'module change' , 'Module Label Span', 'step', 'wn', 'wt'};
    
seqFracListOnlyFracPos = [tempHeader; num2cell(tiedrank(-seqFracListOnlyFracPos(idxTemp,5))), ...
    model.metaboliteName(seqFracListOnlyFracPos(idxTemp,1)),...
    model.metaboliteId(seqFracListOnlyFracPos(idxTemp,1)), ...
    num2cell(seqFracListOnlyFracPos(idxTemp,:))];
%% make rp cell table
mns = mnsResults.mns{1};
for i = 1:size(breakingReaction,1)
    idxRp = find(model.mat(:,mns.results.reactionPairs(i,1)) ~= 0 & model.mat(:,mns.results.reactionPairs(i,2)) ~= 0);
    
    ECString = '';
    geneSymbolString = '';
    rpString = '';
    
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
    rpStruct{i,1} = rpString;
    rpStruct{i,2} = ECString;
    rpStruct{i,3} = geneSymbolString;
    rpStruct{i,4} = [model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}];
    rpStruct{i,5} = idxRp;
end
%% make neighbour fracture outlist
% [~,idxTemp] = sort(nFracList(:,4),'descend');
[~,idxTemp] = sortrows(-nFracList(:,[4 2]));
tempHeader = {'Rank' 'RP ID' 'EC' 'Gene Symbol' 'Reaction' 'Model Rp Id' 'MNS id' '#Fractures'...
    '#Fracture groups' , 'step', 'wn', 'wt'};
tempSum = -nFracList(idxTemp,4)*100-nFracList(idxTemp,2);
nFracList = [tempHeader; num2cell(tiedrank(tempSum)), rpStruct(nFracList(idxTemp,1),:), num2cell(nFracList(idxTemp,:))];

%% make neighbour fracture outlist
[~,idxTemp] = sortrows(-nFracListSplit(:,[5 4]));
tempHeader = {'Rank' 'RP ID' 'EC' 'Gene Symbol' 'Reaction' 'Model Rp Id' 'MNS id' 'StartFrame' 'EndFrame' '#Fractures'...
    , 'step', 'wn', 'wt'};
tempSum = -nFracListSplit(idxTemp,5)*100-nFracListSplit(idxTemp,4);
nFracListSplit = [tempHeader; num2cell(tiedrank(tempSum)), rpStruct(nFracListSplit(idxTemp,1),:), num2cell(nFracListSplit(idxTemp,:))];

end

