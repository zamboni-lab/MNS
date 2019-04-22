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
function met = mns_diffData2ObsEdgeFile(metIdIA, dataStruct, groupMean, ...
    groupStd, filename, initMNS, l1, l2, p2file, nCliques, initLabels)

data = dataStruct.data;
dataAnnotation = dataStruct.annotation;
dataType = dataStruct.dataType;
obsFuncType = initMNS.normObsProb;
initLabelsType = initMNS.initLabelsType;
initLabelsTypeApplyToAll = initMNS.initLabelsTypeApplyToAll;
% obsFuncType
% 1: %normalized gaussian
% 2: %gaussian
% 3: %normalized and linkage dependendt gaussian
% 4: %linkage dependendt gaussian

% link metIdIA to Data ID
aId = [];
aMetId = [];
aMetScore = [];
aMetName = [];
aMetFormula = [];
aNHits = [];
aLink2Ia = []; %first column ID of metIdIA, second column Id of Data (aId)
dataId = [];
c = 0;

% link data annotation 2 model ids
switch dataType
    case 'fiaExp'
        for i = 1:length(dataAnnotation)
            a = dataAnnotation(i);
            for j = 1:length(a.name)
                c = c+1;
                aId = [aId; c];
                aMetId = [aMetId; a.id(j)];
                aMetScore = [aMetScore; a.score(j)];
                aMetName = [aMetName; a.name(j)];
                aMetFormula = [aMetFormula; a.formula(j)];
                aNHits = [aNHits; a.nHits];
                idx = find(strcmp(metIdIA, a.id{j}) == 1);
                dataId = [dataId, i];
                if ~isempty(idx)
                    aLink2Ia = [aLink2Ia; idx c];
                end
            end
        end
    case 'MetIdList'
        for i = 1:length(dataAnnotation)
                c = c+1;
                aId = [aId; i];
                aMetId = [aMetId; dataAnnotation(i)];
                idx = find(strcmp(metIdIA, dataAnnotation{i}) == 1);
                dataId = [dataId, i];
                if ~isempty(idx)
                    aLink2Ia = [aLink2Ia; idx c];
                end
        end
    case 'fiaExp v3.0'
        anndata = dataAnnotation.anndata;
        dataAnnotation = dataAnnotation.annotation;
        for i = 1:length(dataAnnotation)
            a = dataAnnotation(i);
            for j = 1:length(a.anndataIdx)
                c = c+1;
                aId = [aId; c];
                aMetId = [aMetId; anndata.id{a.anndataIdx(j)}];
                aMetScore = [aMetScore; a.score(j)];
                aMetName = [aMetName; anndata.name(a.anndataIdx(j))];
                aMetFormula = [aMetFormula; anndata.formula(a.anndataIdx(j))];
                aNHits = [aNHits; a.nHits];
                idx = find(ismember(metIdIA,anndata.id{a.anndataIdx(j)}));
                dataId = [dataId, i];
                if ~isempty(idx)
                    aLink2Ia = [aLink2Ia; idx c];
                end
            end
        end
    otherwise
        %not defined yet
end

% create vectors
met.modelVarId = zeros(size(metIdIA,1),1);
met.metId = metIdIA;
met.aLink2Ia = aLink2Ia; %links ion
met.detected = zeros(size(metIdIA,1),1);
met.data = zeros(size(metIdIA,1),1);
met.obsCliquePotential = zeros(size(metIdIA,1),length(groupMean)); %not normalized probabilities for each metabolite to be in certain group
met.iaId2dataId = zeros(size(metIdIA,1),1); %links model id 2 data id
met.obsFacId = zeros(size(metIdIA,1),1);

met.dataIonIdx = dataId; %id of ion in fiaExp;
switch initLabelsType
    case 'middle cluster'
        if mod(length(groupMean),2) == 1
            metLabels = ones(size(metIdIA,1),1)*(length(groupMean)-1)/2;
        else
            midLabel = (length(groupMean)-1)/2;
            metLabels = randi([floor(midLabel) ceil(midLabel)], size(metIdIA,1),1);
        end
    case 'infered labels'
        if length(initLabels) == size(metIdIA,1)
            metLabels = initLabels;
        else
            metLabels = randi([0 length(groupMean)-1], size(metIdIA,1),1);
        end
    case 'zeros'
        metLabels = zeros(size(metIdIA,1),1);
    case 'ones'
        metLabels = ones(size(metIdIA,1),1);
    case 'random'
        metLabels = randi([0 length(groupMean)-1], size(metIdIA,1),1);
end

met.initMetLabel = metLabels;
% open outfile for observable edge data
if p2file
    fid = fopen(filename, 'w');
    fprintf(fid,'%d\t%d\r\n', length(metIdIA), length(groupMean));
end
c = 0;
for i = 1:length(metIdIA)
    met.modelVarId(i) = i-1;
    idx = find(aLink2Ia(:,1) == i);
    if ~isempty(idx)
        c = c+1;
        met.obsFacId(i) = c;
        met.detected(i) = 1;
        idTemp = dataId(aLink2Ia(idx,2));
        met.iaId2dataId(i) = idTemp(1); %shoul be as well checked with score
        dataTemp = data(dataId(aLink2Ia(idx,2)));
        if ~initLabelsTypeApplyToAll
            met.initMetLabel(i) = initLabels(dataId(aLink2Ia(idx,2)));
        end
        met.data(i) = dataTemp(1); %Should be compared with score
        met.obsCliquePotential(i,:) = mns_calcProbGaussian(groupMean, groupStd,...
            met.data(i),obsFuncType, l2, nCliques.nCliquePerVar(i,:),nCliques.nNeighboorMetabolites(i,:),l1);
    end
    if p2file
        for j = -2:length(groupMean)
            if j == -2
                fprintf(fid,'%d', met.modelVarId(i));
            elseif j == -1
                fprintf(fid,'\t%d', met.detected(i));
            elseif j == 0
                fprintf(fid,'\t%d', met.initMetLabel(i));
            else
                fprintf(fid,'\t%f', met.obsCliquePotential(i,j));
            end
            
        end
        
        if i ~= length(metIdIA)
            fprintf(fid,'\r\n');
        end
    end
end

if p2file
    fclose(fid);
end


end

% % calc probability based on gaussian
% 
% function [prob] = calcProbGaussian(meanVec, stdVec, data, norm, l2)
% prob = zeros(length(meanVec),1);
% 
% for i = 1:length(prob)
%     if norm
%         prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( - l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%     else
%         prob(i) = exp( - l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%     end
% end
% 
% end