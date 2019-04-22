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
function mns_plotTcWithModules(mnsResults, mnsDataStruct, model, idx, plotMode, metIdx)
if nargin < 5
    plotMode = 'area';
end

tempResults = mnsResults.mns{1};
nTp = size(mnsDataStruct.data,2);
if nargin < 6
    metIdx = 1:length(mnsResults.mns{1}.observations.detected)/nTp;
end
nMetData = size(mnsDataStruct.data,1);
nMetModel = length(metIdx);
% uTime = unqiue
clusterLabel = mnsResults.clustLabel(:,:,idx);
% cmap = lines(max(max(clusterLabel))+1);
% cmap = redbluecmap(max(max(clusterLabel))+1);
%%
maxCl = max(max(clusterLabel));
cSize = maxCl*50;
cmap = [];
cmapTemp = [];

switch plotMode
    case 'area'
        cwhite= [247 247 247]/255;
    case 'line'
        cwhite= [180 180 180]/255;
end

cRef = [0 102 204]/255*0.8;
cpos = cwhite;

cDelta = cpos-cRef;
for i = 1:3
    if cDelta(i) ~= 0
        cmapTemp(1:cSize,i) = cRef(i):cDelta(i)/(cSize-1):cpos(i);
    else
        cmapTemp(1:cSize,i) = repmat(cpos(i),cSize,1);
    end
end
cmap = [cmap;cmapTemp];
cpos = [204 0 51]/255;
cRef = cwhite;
cDelta = cpos-cRef;
for i = 1:3
    if cDelta(i) ~= 0
        cmapTemp(1:cSize,i) = cRef(i):cDelta(i)/(cSize-1):cpos(i);
    else
        cmapTemp(1:cSize,i) = repmat(cpos(i),cSize,1);
    end
end
cmap = [cmap;cmapTemp];
if maxCl <= 2
    cmap = cmap([50:50:end-50 end-50],:);
else
    cmap = cmap([1:100:end end],:);
end
%%
if isfield(mnsDataStruct, 'time')
    tp = mnsDataStruct.time;
else
    tp = 1:nTp;
end

%% sum up the data
dataTemp = zeros(nMetData, nTp+1);
for i = 2:nTp+1
    dataTemp(:,i) = mnsDataStruct.data(:,i-1)+dataTemp(:,i-1);
end
%%

nW = 5;
nH = 4;
if nMetModel < nW*nH
    nH = floor(sqrt(nMetModel));
    if nMetModel > nW*nH
        nW = nH;
    else
        nW = nH+1;
    end
end
%%
c = 0;
figure()
for i = 1:nMetModel
    c = c+1;
    if c > nW * nH
        c = 1;
        figure();
    end
    subplot(nH, nW, c);
%     idObservation = (i-1)*nTp+1;
    idObservation = (metIdx(i)-1)*nTp+1;
    ymax = 1;
    ymin = 0;
    idData = [];
    if tempResults.observations.detected(idObservation)
        idData = tempResults.observations.iaId2dataId(idObservation);
%         disp(idData)
        if min(dataTemp(idData,:)) < ymin
           ymin = 1.1* min(dataTemp(idData,:));
        end
        if max(dataTemp(idData,:)) > ymax
           ymax = 1.1* max(dataTemp(idData,:));
        end
    end
    %% plot the Data
    dy = ymax-ymin;
    hold all;
    switch plotMode
        case 'area'
            for j = 1:(length(tp))
%                 clLabel = clusterLabel(i,j)+1;
                clLabel = clusterLabel(metIdx(i),j)+1;
                rectangle('Position',[j ymin 1 dy], 'FaceColor',cmap(clLabel,:),'EdgeColor',[0.3 0.3 0.3], 'LineWidth', 0.5)
%                 if ~isempty(idData)
%                     %                 plot(tp, mnsDataStruct.data(idData,:), '-k');
%                     plot([j j+1], dataTemp(idData,[j j+1]), '-k', 'LineWidth', 2);
%                 end
            end
            if ~isempty(idData)
                %                 plot(tp, mnsDataStruct.data(idData,:), '-k');
                plot(1:length(tp)+1, dataTemp(idData,:), '-k', 'LineWidth', 2);
            end
            
        case 'line'
            if ~isempty(idData)
                for j = 1:(length(tp))
%                     clLabel = clusterLabel(i,j)+1;
                    clLabel = clusterLabel(metIdx(i),j)+1;
                    hold all
                    plot([j j+1], dataTemp(idData,[j j+1]), '-', 'Color', cmap(clLabel,:), 'LineWidth', 2);
                end
            else
                for j = 1:(length(tp))
%                     clLabel = clusterLabel(i,j)+1;
                    clLabel = clusterLabel(metIdx(i),j)+1;
%                     hold all
                    plot([j j+1], [0.5 0.5], '-', 'Color', cmap(clLabel,:), 'LineWidth', 2);
                end
            end
    end
    
    xlabel('timepoint');
    ylabel('log2fc');
    title(model.metaboliteName(metIdx(i)));
    xlim([min(tp) max(tp)+1]);
    ylim([ymin ymax]);
end


end