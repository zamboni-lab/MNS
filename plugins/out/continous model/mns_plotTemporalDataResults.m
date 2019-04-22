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
function mns_plotTemporalDataResults(mnsResults,model,idx, figLabel)
if nargin < 4
    label = '';
end
% idx = find(penArr == max(max(penArr)));
% idx = idx(1);
% maxTL3 = ceil(idx/size(penArr,1));
% maxNl1 = mod(idx, size(penArr,1));
% if maxNl1 == 0
%     maxNl1 = size(penArr,1);
% end
% idx = find(mnsResults.nL1 == uNL1(maxNl1) & mnsResults.tL3 == uTL3(maxTL3));
figure('Name', figLabel')
subplot(1,3,1)
hold all
% tempidx = mnsResults.mns{idx}.observations.modelVarIdTimePoint+1;
tempidx = mnsResults.mns{1}.observations.modelVarIdTimePoint+1;
% switch max(max(mnsResults.mns{idx}.results.inferedLabels(tempidx)))+1
%     case 2
%         cmap = [130 204 244;244 122 215]/255;
%     case 3
% %         cmap = [0 174 239;247 247 247;248 46 192]/255;
%         cmap = [130 204 244;247 247 247;244 122 215]/255;
%     case 4
%         cmap = [0 174 239;130 204 244;244 122 215;248 46 192]/255;
%     case 5
%         cmap = [0 174 239;130 204 244;247 247 247;244 122 215;248 46 192]/255;
%     otherwise
%         cmap = redbluecmap(max(max(max(max(mnsResults.mns{idx}.results.inferedLabels(tempidx)))+1))+1);
% end
%%
% maxCl = max(max(mnsResults.mns{idx}.results.inferedLabels(tempidx)));
maxCl = max(max(mnsResults.clustLabel(:,:,idx)));
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

%%
% cmap = redbluecmap(max(max(mnsResults.mns{idx}.results.inferedLabels(tempidx)))+1);
imagesc(mnsResults.clustLabel(:,:,idx))
% imagesc(mnsResults.mns{idx}.results.inferedLabels(tempidx))
set(gca,'YDir','normal')
set(gca, 'YTick', 1:length(model.metaboliteName), 'YTickLabel', model.metaboliteName)
xlim([0.5 size(tempidx,2)+0.5])
ylim([0.5 size(tempidx,1)+0.5])
title('Module labels');
xlabel('frame');
set(gca, 'CLim', [0 mnsResults.mns{1}.parameters.noOfLabels-1])
subplot(1,3,2)
imagesc(mnsResults.temporalBreaks(:,:,idx)+1, [0 2])
title(sprintf('Sequential fractures; n=%d',sum(sum(mnsResults.temporalBreaks(:,:,idx)))));
set(gca, 'XTick', 1:1:size(tempidx,2)-1, 'XTickLabel', 1.5:1:size(tempidx,2)-0.5)
set(gca, 'YTick', 1:length(model.metaboliteName), 'YTickLabel', model.metaboliteName)
% xlim([0.5 size(tempidx,2)+0.5])
xlabel('frame intersect');
set(gca,'YDir','normal')

%%
mns = mnsResults.mns{1};
for i = 1:length(mns.results.reactionPairs)
    if mns.results.reactionPairs(i,1) > length(model.metaboliteName)
        break;
    end
    rLabel(i) = {[model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}]};
end
%%
subplot(1,3,3)
imagesc(mnsResults.breakingReaction(:,:,idx)+1, [0 2])
title(sprintf('Neighborhood fractures; n=%d',sum(sum(mnsResults.breakingReaction(:,:,idx)))));
set(gca, 'YTick', 1:length(rLabel), 'YTickLabel', rLabel)
xlabel('frame')
set(gca,'YDir','normal')
colormap(cmap)


end