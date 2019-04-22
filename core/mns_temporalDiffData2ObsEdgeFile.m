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
function met = mns_temporalDiffData2ObsEdgeFile(metIdIA, data, dataAnnotation, dataType, ...
    groupMean, groupStd, filename, initMNS, l2, p2file, nCliques, timeFileName, tL3, initLabels)

% obsFuncType
% 1: %normalized gaussian
% 2: %gaussian
% 3: %normalized and linkage dependendt gaussian
% 4: %linkage dependendt gaussian

obsFuncType = initMNS.normObsProb;
initLabelsType = initMNS.initLabelsType;
initLabelsTypeApplyToAll = initMNS.initLabelsTypeApplyToAll;


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
nTp = size(data,2);
offSet = size(nCliques.varId(:,1),1);

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
                idx = find(strcmp(metIdIA, anndata.id{a.anndataIdx(j)}{:}) == 1);
                dataId = [dataId, i];
                if ~isempty(idx)
                    aLink2Ia = [aLink2Ia; idx c];
                end
            end
        end
    otherwise
end

% create vectors
met.modelVarId = zeros(size(metIdIA,1)*nTp,1);
met.metId = metIdIA;
met.detected = zeros(size(metIdIA,1)*nTp,1);
met.data = zeros(size(metIdIA,1)*nTp,1);
met.aLink2Ia = aLink2Ia;
met.obsCliquePotential = zeros(size(metIdIA,1)*nTp,size(groupMean,1)); %not normalized probabilities for each metabolite to be in certain group
met.iaId2dataId = zeros(size(metIdIA,1)*nTp,1);
met.obsFacId = zeros(size(metIdIA,1)*nTp,1);
met.noOfTimePoints = nTp;

met.dataIonIdx = dataId; %id of ion in fiaExp;
switch initLabelsType
    case 'middle cluster'
        if mod(length(groupMean),2) == 1
            metLabels = ones(size(metIdIA,1)*nTp,1)*(length(groupMean)-1)/2;
        else
            midLabel = (size(groupMean,1)-1)/2;
            metLabels = randi([floor(midLabel) ceil(midLabel)], size(metIdIA,1)*nTp,1);
        end
    case 'infered labels'
        if length(initLabels) == size(metIdIA,1)*nTp
            metLabels = initLabels;
        else
            metLabels = randi([0 length(groupMean)-1], size(metIdIA,1)*nTp,1);
        end
    case 'zeros'
        metLabels = zeros(size(metIdIA,1)*nTp,1);
    case 'ones'
        metLabels = ones(size(metIdIA,1)*nTp,1);
    case 'random'
        metLabels = randi([0 length(groupMean)-1], size(metIdIA,1)*nTp,1);
end

met.initMetLabel = metLabels;

c = 0;
for i = 1:length(metIdIA)
    for j = 1:nTp
        c = c+1;
        met.metId{c} = metIdIA{i};
    end
end

% open outfile for observable edge data
if p2file
    fid = fopen(filename, 'w');
    fidTime = fopen(timeFileName, 'w');
    fprintf(fid,'%d\t%d\r\n', length(metIdIA)*size(data,2), size(groupMean,1));
    %timeFile: number of rows; number of time windows, potential if not same
    %label, potential if same label;
    fprintf(fidTime,'%d\t%d\t%f\t%f\r\n', length(metIdIA)*(nTp-1), nTp, exp(-tL3), exp(tL3));
end
c = 0;

for i = 1:length(metIdIA)
    for k = 1:nTp
        id = (i-1)*nTp+k;
%         if k > 1
%             fprintf(fidTime,'%d\t', (k-2)*offSet+i-1);
%         end
        met.modelVarId(id) = (k-1)*offSet+i-1;
        if k > 1
            fprintf(fidTime,'%d\t%d', (k-2)*offSet+i-1, met.modelVarId(id));
        end
        
        met.modelVarIdTimePoint(i,k) = (k-1)*offSet+i-1;
        idx = find(aLink2Ia(:,1) == i);
%         fprintf(fidTime,'%d\t%d\r\n', met.modelVarId(id));
        if ~isempty(idx)
            c = c+1;
            met.obsFacId(id) = c;

            met.detected(id) = 1;
            dataTemp = data(dataId(aLink2Ia(idx,2)),k);
            if ~(initLabelsTypeApplyToAll)
                met.initMetLabel(id) = initLabels(dataId(aLink2Ia(idx,2)),k);
            end
            met.iaId2dataId(id) = dataId(aLink2Ia(idx(1),2));
            met.data(id) = dataTemp(1); %Should be compared with score
            met.obsCliquePotential(id,:) = mns_calcProbGaussian(groupMean(:,k), groupStd(:,k), met.data(id),...
                obsFuncType, l2, nCliques.nCliquePerVar(:,k),nCliques.nNeighboorMetabolites(i,:));
        end
        if p2file
            for j = -2:size(groupMean,1)
                if j == -2
                    fprintf(fid,'%d', met.modelVarId(id));
                elseif j == -1
                    fprintf(fid,'\t%d', met.detected(id));
                elseif j == 0
                    fprintf(fid,'\t%d', met.initMetLabel(id));
                else
                    fprintf(fid,'\t%f', met.obsCliquePotential(id,j));
                end
                
            end
            if ~(i == length(metIdIA) && k == nTp)
                fprintf(fid,'\r\n');
                if k > 1
                    fprintf(fidTime,'\r\n');
                end
            end
        end
    end
end
if p2file
    fclose(fid);
    fclose(fidTime);
end


end

% % calc probability based on gaussian
% 
% function [prob] = calcProbGaussian(meanVec, stdVec, data, obsFuncType,l2)
% prob = zeros(length(meanVec),1);
% 
% for i = 1:length(prob)
%     switch obsFuncType
%         case 1 %normalized gaussian
%             prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%         case 2 %gaussian
%             prob(i) = exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%         case 3 %normalized and linkage dependendt gaussian
%             prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%         case 4 %linkage dependendt gaussian
%             prob(i) = exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%     end
% end
% 
% end