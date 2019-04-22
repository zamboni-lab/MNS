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
function met = mns_diffDataBiolRep2ObsEdgeFile(metIdIA, data, dataAnnotation, dataType,...
    groupMean, groupStd, filename, obsFuncType, l2, p2file, nCliques, avgType)

if nargin < 12 || isempty(avgType)
    avgType = 'median';
end
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
nRepl = size(data,2);
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
    otherwise
end
% create vectors
met.modelVarId = zeros(size(metIdIA,1)*nRepl,1);
met.metId = metIdIA;
met.detected = zeros(size(metIdIA,1)*nRepl,1);
met.data = zeros(size(metIdIA,1)*nRepl,1);
met.aLink2Ia = aLink2Ia;
met.obsCliquePotential = zeros(size(metIdIA,1)*nRepl,length(groupMean)); %not normalized probabilities for each metabolite to be in certain group
met.iaId2dataId = zeros(size(metIdIA,1),1);
met.obsFacId = zeros(size(metIdIA,1),1);
c = 0;
for i = 1:length(metIdIA)
    for j = 1:nRepl
        c = c+1;
        met.metId(c) = metIdIA(i);

    end
end

% open outfile for observable edge data
if p2file
    fid = fopen(filename, 'w');
    fprintf(fid,'%d\t%d\r\n', length(metIdIA)*size(data,2), length(groupMean));
end
c = 0;

% old way
% for i = 1:length(metIdIA)
%     for k = 1:nRepl
%         id = (i-1)*nRepl+k;
%         met.modelVarId(id) = i-1;
%         idx = find(aLink2Ia(:,1) == i);
%         
%         if ~isempty(idx)
%             c = c+1;
%             met.obsFacId(id) = c;
% 
%             met.detected(id) = 1;
%             dataTemp = data(dataId(aLink2Ia(idx,2)),k);
%             met.iaId2dataId(id) = dataId(aLink2Ia(idx,2));
%             met.data(id) = dataTemp(1); %Should be compared with score
%             met.obsCliquePotential(id,:) = mns_calcProbGaussian(groupMean, groupStd, met.data(id), ...
%                 obsFuncType, l2,nCliques.nCliquePerVar);
%         end
%         if p2file
%             for j = -1:length(groupMean)
%                 if j == -1
%                     fprintf(fid,'%d', met.modelVarId(id));
%                 elseif j == 0
%                     fprintf(fid,'\t%d', met.detected(id));
%                 else
%                     fprintf(fid,'\t%f', met.obsCliquePotential(id,j));
%                 end
%                 
%             end
%             if i ~= length(metIdIA) || k ~= nRepl
%                 fprintf(fid,'\r\n');
%             end
%         end
%     end
% end

%NEW WAY
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
        obsCliquePotential = zeros(nRepl,length(groupMean));
        dataTempArr = zeros(nRepl,1);
        for k = 1:nRepl
            dataTemp = data(dataId(aLink2Ia(idx,2)),k);
            dataTempArr(k) = dataTemp(1,:); %Should be compared with score
            obsCliquePotential(k,:) = mns_calcProbGaussian(groupMean, groupStd,dataTempArr(k),...
                obsFuncType, l2, nCliques.nCliquePerVar);
        end
        switch avgType
            case 'mean'
                met.obsCliquePotential(i,:) = mean(obsCliquePotential,1);
                met.data(i) = mean(dataTempArr(k));
            case 'median'
                met.obsCliquePotential(i,:) = median(obsCliquePotential,1);
                met.data(i) = median(dataTempArr(k));
        end
        clear dataTempArr obsCliquePotential;
    end
    if p2file
        for j = -1:length(groupMean)
            if j == -1
                fprintf(fid,'%d', met.modelVarId(i));
            elseif j == 0
                fprintf(fid,'\t%d', met.detected(i));
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
% function [prob] = calcProbGaussian(meanVec, stdVec, data, norm,l2)
% prob = zeros(length(meanVec),1);
% 
% for i = 1:length(prob)
%     if norm
%         prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%     else
%         prob(i) = exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
%     end
% end
% 
% end