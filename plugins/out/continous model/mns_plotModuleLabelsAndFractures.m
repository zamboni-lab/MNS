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
function mns_plotModuleLabelsAndFractures(mnsResults, model,wsVec,wnVec)
%% 
nSubplot = length(wsVec);
%% make cmap
maxCl = max(max(max(mnsResults.clustLabel(:,:,:))));
cSize = maxCl*50;
cmap = [];
cmapTemp = [];

cRef = [0 102 204]/255*0.8;
cpos= [247 247 247]/255;
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
cRef = [247 247 247]/255;
cDelta = cpos-cRef;
for i = 1:3
    if cDelta(i) ~= 0
        cmapTemp(1:cSize,i) = cRef(i):cDelta(i)/(cSize-1):cpos(i);
    else
        cmapTemp(1:cSize,i) = repmat(cpos(i),cSize,1);
    end
end
cmap = [cmap;cmapTemp];
cmap = cmap([1:100:end end],:);

%% get idx and plot module labels
figure('Name','Module labels')
for i = 1:length(wsVec)
    [~,~, idx(i)] = mns_calcProbabilityScanData(mnsResults, model,wsVec(i),wnVec(i),false, false);
    subplot(1,nSubplot,i)
    hold all
    tempidx = mnsResults.mns{1}.observations.modelVarIdTimePoint+1;
    
        %%
    % cmap = redbluecmap(max(max(mnsResults.mns{idx}.results.inferedLabels(tempidx)))+1);
    imagesc(mnsResults.clustLabel(:,:,idx(i)))
    % imagesc(mnsResults.mns{idx}.results.inferedLabels(tempidx))
    set(gca,'YDir','normal')
    if i == 1
        set(gca, 'YTick', 1:length(model.metaboliteName), 'YTickLabel', model.metaboliteName)
    else
        set(gca, 'YTick', 1:length(model.metaboliteName), 'YTickLabel', '')
    end
    xlim([0.5 size(tempidx,2)+0.5])
    ylim([0.5 size(tempidx,1)+0.5])
    title(sprintf('ws = %0.2f, wn = %0.2f', wsVec(i), wnVec(i)));
    xlabel('frame');
    set(gca, 'CLim', [0 mnsResults.mns{1}.parameters.noOfLabels-1])
end
colormap(cmap)
%% plot sequential fractures
figure('Name','Sequential fractures')
for i = 1:length(idx)
    subplot(1,nSubplot,i)
    imagesc(mnsResults.temporalBreaks(:,:,idx(i))+1, [0 2])
    title(sprintf('ws = %0.2f, wn = %0.2f; n=%d', wsVec(i), wnVec(i),sum(sum(mnsResults.temporalBreaks(:,:,idx(i))))));
    set(gca, 'XTick', 1:1:size(tempidx,2)-1, 'XTickLabel', 1.5:1:size(tempidx,2)-0.5)
    if i == 1
        set(gca, 'YTick', 1:length(model.metaboliteName), 'YTickLabel', model.metaboliteName)
    else
        set(gca, 'YTick', 1:length(model.metaboliteName), 'YTickLabel', '')   
    end
    % xlim([0.5 size(tempidx,2)+0.5])
    xlabel('frame intersect');
    set(gca,'YDir','normal')
end
colormap(cmap)
%% plot neighborhood fractures
mns = mnsResults.mns{1};
for i = 1:length(mns.results.reactionPairs)
    if mns.results.reactionPairs(i,1) > length(model.metaboliteName)
        break;
    end
    rLabel(i) = {[model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}]};
end

%
figure('Name','Neighborhood fractures')
for i = 1:length(idx)
    subplot(1,nSubplot,i)
    imagesc(mnsResults.breakingReaction(:,:,idx(i))+1, [0 2])
    title(sprintf('ws = %0.2f, wn = %0.2f;  n=%d', wsVec(i), wnVec(i),sum(sum(mnsResults.breakingReaction(:,:,idx(i))))));
    if i == 1
        set(gca, 'YTick', 1:length(rLabel), 'YTickLabel', rLabel)
    else
        set(gca, 'YTick', 1:length(rLabel), 'YTickLabel', '')
    end
    xlabel('frame')
    set(gca,'YDir','normal')

end
colormap(cmap)
end
