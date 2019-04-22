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
function mnsResultsSummary = mns_plotBreakingReactions(mnsResults, metName)
mns = mnsResults{1};
mode = mns.results.mode;

switch mode
    case 'two-state'
        mnsResultsSummary = twoStateFractures(mnsResults, metName);
    case 'timeframe'
        
end


end

function mnsResultsSummary = temporalModelFractures(mnsResults, metName)


end


% calculate, plot and summarize fracture information of a two state MRF
% Model
function mnsResultsSummary = twoStateFractures(mnsResults, metName)

mns = mnsResults{1};
for i = 1:length(mnsResults)
    breakingReaction(:,i) = mnsResults{i}.results.breakingReactions;
    clustLabel(:,i) = mnsResults{i}.results.inferedLabels;
    nL(i) = mnsResults{i}.parameters.nL1;
end

totalFracturesPerReaction = sum(breakingReaction, 2);
totalFracturesPerReactionRank = tiedrank(-totalFracturesPerReaction);
[~,~,totalFracturesPerReactionRank] = unique(totalFracturesPerReactionRank);

idVec = 1:size(breakingReaction,2);
fractureIdPerReaction = zeros(size(breakingReaction));

for i = 1:size(breakingReaction,1)
    fractureIdPerReaction(i,breakingReaction(i,:) == 1) = idVec(breakingReaction(i,:) == 1);
end

% max nL1 ID for each RP where a fracture is occuring
maxFractureIdForRP = max(fractureIdPerReaction, [], 2);
maxFractureIdForRpRank = tiedrank(-maxFractureIdForRP);
[~,~,maxFractureIdForRpRank] = unique(maxFractureIdForRpRank);

figure, imagesc(breakingReaction)
for i = 1:size(mns.results.reactionPairs,1)
    label{i} = [metName{mns.results.reactionPairs(i,1)} ' <-> ' metName{mns.results.reactionPairs(i,2)}];
end

colormap([1 1 1; 0.75 0 0])
title('Fractured reactions')
set(gca, 'YTick', 1:size(mns.results.reactionPairs,1), 'YTickLabel', label);
set(gca, 'XTick', 1:length(nL), 'XTickLabel', nL);
xlabel('\lambda_1')

% plot total fractures
totalBreaks = sum(breakingReaction,1);
figure, plot(nL, totalBreaks, '-', 'LineWidth', 2);
title('Total fractures');
xlabel('\lambda_1');
ylabel('# of fractures');

%plot change in fractures
diffBreaks = zeros(size(totalBreaks));
diffBreaks(2:end) = totalBreaks(1:end-1)-totalBreaks(2:end);
figure, plot(nL, diffBreaks, '-', 'LineWidth', 2);
title('Change in # of fractures');
xlabel('\lambda_1');
ylabel('\Delta # of fractures');

figure, imagesc(clustLabel)

% for i = 1:size(mns.results.reactionPairs,1)
%     label{i} = [model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}];
% end
% colormap([1 1 1; 0.75 0 0])

title('ClusterLabels')
set(gca, 'YTick', 1:size(mns.results.reactionPairs,1), 'YTickLabel', metName);
set(gca, 'XTick', 1:length(nL), 'XTickLabel', nL);
xlabel('\lambda_1')
cmap = colormap(jet(mnsResults{1}.parameters.noOfLabels));
cmap = [236 0 140; 46 49 146; 0 174 239; 0 166 81;255 242 0]/255;
colormap(cmap)

mnsResultsSummary.metabolite = metName;
mnsResultsSummary.reaction = label;
mnsResultsSummary.clusterLabels = clustLabel;
mnsResultsSummary.fractureDetailed = breakingReaction;
mnsResultsSummary.totalBreaks = totalBreaks;
mnsResultsSummary.diffBreaks = diffBreaks;
mnsResultsSummary.totalFracturesPerReaction = totalFracturesPerReaction;
mnsResultsSummary.totalFracturesPerReactionRank = totalFracturesPerReactionRank;
mnsResultsSummary.maxFractureIdForRP = maxFractureIdForRP;
mnsResultsSummary.maxFractureIdForRpRank = maxFractureIdForRpRank;

end