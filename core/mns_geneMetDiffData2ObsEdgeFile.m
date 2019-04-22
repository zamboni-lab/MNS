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
function out = mns_geneMetdiffData2ObsEdgeFile(dataStruct, model, groupMean, ...
    groupStd, filename, initMNS, l1, l2, nCliques, initLabels, metObservations, p2file)
if nargin < 11
    metObservations = [];
end
if nargin < 12
    p2file = true;
end
data = dataStruct.geneData;
dataAnnotation = dataStruct.geneAnnotation;
dataType = dataStruct.geneDataType;
obsFuncType = initMNS.normObsProb;
initLabelsType = initMNS.initLabelsType;
initLabelsTypeApplyToAll = initMNS.initLabelsTypeApplyToAll;
% obsFuncType
% 1: %normalized gaussian
% 2: %gaussian
% 3: %normalized and linkage dependendt gaussian
% 4: %linkage dependendt gaussian

%% link geneData to model gene and Rp ID

% aLink2Ia = []; %first column ID of metIdIA, second column Id of Data (aId)
dataId = [];
a2modelGeneId = [];
dataIdx2modelRpId = [];
a2modelRpId = [];

%% link data annotation 2 model ids
switch dataType
    case 'GeneSymbol'
        % get rpId 2 geneId mapping
        for i = 1:length(dataAnnotation)
            
%             a = dataAnnotation(i);
%             geneAnnotationId = [geneAnnotationId; c];
            idx = find(strcmp(model.GeneSymbol, dataAnnotation{i}) == 1);
            a2modelGeneId = [a2modelGeneId; idx];
            if ~isempty(idx)
                idxRp = find(model.gene2rp(idx,:) == 1);
                a2modelRpId = [a2modelRpId; {idxRp}];
                for j = 1:length(idxRp)
                    dataIdx2modelRpId = [dataIdx2modelRpId; i idxRp(j)];
                end
            end

        end
    otherwise
        disp('not yet implemented');
end


%% generate initial labels
switch initLabelsType
    case 'middle cluster'
        if mod(length(groupMean),2) == 1
            geneLabels = ones(length(model.rpId),1)*(length(groupMean)-1)/2;
        else
            midLabel = (length(groupMean)-1)/2;
            geneLabels = randi([floor(midLabel) ceil(midLabel)], length(model.rpId),1);
        end
    case 'infered labels'
        if length(initLabels) == size(metIdIA,1)
            geneLabels = initLabels;
        else
            geneLabels = randi([0 length(groupMean)-1], length(model.rpId),1);
        end
    case 'zeros'
        geneLabels = zeros(length(model.rpId),1);
    case 'ones'
        geneLabels = ones(length(model.rpId),1);
    case 'random'
        geneLabels = randi([0 length(groupMean)-1], length(model.rpId),1);
end
%% prepare out structure
noOfRp = length(model.rpId);
noOfVariables = noOfRp;
out.initLabel = geneLabels;
out.isRp = model.isRp;
out.isMet = model.isMet;
if ~isempty(metObservations)
   noOfVariables = noOfVariables+length(metObservations.modelVarId);
   out.initLabel = [out.initLabel;metObservations.initMetLabel];
   out.met = metObservations;
   out.met.aLink2Ia(:,1) = out.met.aLink2Ia(:,1)+noOfRp;
   out.met.modelVarId = out.met.modelVarId+noOfRp;
end
out.modelVarId = zeros(noOfVariables,1);
dataId = 1:length(dataStruct.geneData);
out.noOfVariables = noOfVariables;
out.noOfRp = noOfRp;
% out.aLink2Ia = aLink2Ia; %links ion
out.detected = zeros(noOfVariables,1);
out.data = zeros(noOfVariables,1);
out.obsCliquePotential = zeros(noOfVariables,length(groupMean)); %not normalized probabilities for each metabolite to be in certain group
out.iaId2dataId = cell(noOfVariables,1); %links model id 2 data id
out.obsFacId = zeros(noOfVariables,1);
%% open outfile for observable edge data
if p2file
    fid = fopen(filename, 'w');
    fprintf(fid,'%d\t%d\r\n', noOfVariables, length(groupMean));
end
c = 0;
for i = 1:noOfVariables
    out.modelVarId(i) = i-1;
    if out.isRp(i) && ~out.isMet(i)
        idx = find(dataIdx2modelRpId(:,2) == i);
        if ~isempty(idx)
            c = c+1;
            out.obsFacId(i) = c;
            out.detected(i) = 1;
            idTemp = dataId(dataIdx2modelRpId(idx,1));
            dataTemp = data(dataIdx2modelRpId(idx,1));
            % this so far works only for the best option... and not all
            % options
            [out.data(i), out.iaId2dataId(i), out.dataLabel(i)] = parseGeneData(idTemp,dataTemp,dataAnnotation(idTemp),initMNS);
            if ~initLabelsTypeApplyToAll
                out.initLabel(i) = initLabels(out.iaId2dataId{i}(1));
            end
            out.obsCliquePotential(i,:) = mns_calcProbGaussian(groupMean, groupStd,...
                out.data(i),obsFuncType, l2, nCliques.nCliquePerVar(i,:),nCliques.nNeighboorMetabolites(i,:),l1);
        end
    elseif ~out.isRp(i) && out.isMet(i)
        c = c+1;
        metIdx = (i-noOfRp);
        out.obsFacId(i) = c;
        out.detected(i) = out.met.detected(metIdx);
        out.initLabel(i) = out.met.initMetLabel(metIdx);
        out.data(i) = out.met.data(metIdx);
        out.iaId2dataId(i) = num2cell(out.met.iaId2dataId(metIdx));
        out.dataLabel(i) = out.met.metId(metIdx);
        out.obsCliquePotential(i,:) = out.met.obsCliquePotential(metIdx,:);
    end
    if p2file
        for j = -2:length(groupMean)
            if j == -2
                fprintf(fid,'%d', out.modelVarId(i));
            elseif j == -1
                fprintf(fid,'\t%d', out.detected(i));
            elseif j == 0
                fprintf(fid,'\t%d', out.initLabel(i));
            else
                fprintf(fid,'\t%f', out.obsCliquePotential(i,j));
            end
            
        end
        
        if i ~= noOfVariables
            fprintf(fid,'\r\n');
        end
    end
end

if p2file
    fclose(fid);
end

% out.dataMetIonIdx = dataId; %id of ion in fiaExp;
end

function [outData, outIdx, outGene] = parseGeneData(idx,data,geneSymbol,initMNS)

uGeneSymbol = unique(geneSymbol);
outGene = uGeneSymbol;
%%
for i = 1:length(unique(geneSymbol))
    idxGene = find(strcmp(geneSymbol, uGeneSymbol{i}) == 1);
    switch initMNS.subGraphExtraction.geneAvgType
        case 'all'
            outGene = geneSymbol;
            outData = data;
            outIdx = num2cell(idx);
            break;
        case 'median'
            outData(i) = median(data(idxGene));
            idxOutTemp = find(data(idxGene) == outData(i));
            outIdx(i) = {idx(idxGene)};
        case 'mean'
            outData(i) = mean(data(idxGene));
            outIdx(i) = {idx(idxGene)};
        case 'max'
            outData(i) = max(abs(data(idxGene)));
            outData(i) = sign(data(abs(data(idxGene))==outData(i)))*outData(i);
            outIdx(i) = {idx(idxGene(data(idxGene) == outData(i)))};
    end
    
end
%%
switch initMNS.subGraphExtraction.geneObservations
    case 'best'
        [outData, idxOutBest] = max(outData);
        outIdx = outIdx(idxOutBest);
        outGene = outGene(idxOutBest);
    case 'all'
        disp('needs to be implemented')
end

end