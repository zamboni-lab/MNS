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
function mns_plotTemporalScreeningResults(mnsResults,probArr, score, plotScoreOnly)

plotProbPen = true;
spSize = 4;
if nargin < 2
    plotProbPen = false;
    spSize = 2;
end
if nargin < 4
    plotScoreOnly = false;
end
scoreTemp = score - nanmin(nanmin(score));
scoreTemp = scoreTemp/nanmax(nanmax(scoreTemp));

if plotScoreOnly
%     imagesc(score./max(max(score)), [0.9 1]);
    imagesc(scoreTemp, [0 1]);
    set(gca,'YDir','normal')
    xlim([0.5 size(score,2)+0.5])
    ylim([0.5 size(score,1)+0.5])
    ylabel('nL1 range');
    xlabel('tL3 range');
    set(gca, 'XTick', 1:10:length(mnsResults.tL3Range), 'XTickLabel', round(mnsResults.tL3Range(1:10:length(mnsResults.tL3Range)),3))
    set(gca, 'YTick', 1:10:length(mnsResults.nL1Range), 'YTickLabel', round(mnsResults.nL1Range(1:10:length(mnsResults.nL1Range)),3))
    % set(gca,'YTick', 1:length(uNL1),'YTickLabel', nL1labels);
    % set(gca,'XTick', 1:length(uTL3),'XTickLabel', uTL3);
    title('score');
    idx = find(score == max(max(score)));
    for i = 1:length(idx)
        idxTemp = idx(i);
        maxTL3 = ceil(idxTemp/size(score,1));
        maxNl1 = mod(idxTemp, size(score,1));
        if maxNl1 == 0
            maxNl1 = size(score,1);
        end
        text(maxTL3-0.2,maxNl1-0.1,'*', 'FontSize', 10)
    end
    colormap(makeCmap);
else
    figure
    subplot(1,spSize,1)
    imagesc(mnsResults.amountBreakingReaction./max(max(mnsResults.amountBreakingReaction)), [0 1]);
    set(gca,'YDir','normal')
    ylabel('\lambda_1');
    xlabel('\lambda_2');
    % set(gca,'YTick', 1:length(uNL1),'YTickLabel', nL1labels);
    % set(gca,'XTick', 1:length(uTL3),'XTickLabel', uTL3);
    title('Sum of neighborhood fractures');
    colorbar
    colormap(makeCmap)
    %
    subplot(1,spSize,2)
    imagesc(mnsResults.amountTemporalBreaks./max(max(mnsResults.amountTemporalBreaks)), [0 1]);
    set(gca,'YDir','normal')
    ylabel('\lambda_1');
    xlabel('\lambda_2');
    % set(gca,'YTick', 1:length(uNL1),'YTickLabel', nL1labels);
    % set(gca,'XTick', 1:length(uTL3),'XTickLabel', nL1labels);
    title('Sum of sequential fractures');
    %
    if plotProbPen
        subplot(1,spSize,3)
        imagesc(probArr./max(max(probArr)), [0.9 1]);
        set(gca,'YDir','normal')
        ylabel('\lambda_1');
        xlabel('\lambda_2');
        % set(gca,'YTick', 1:length(uNL1),'YTickLabel', nL1labels);
        % set(gca,'XTick', 1:length(uTL3),'XTickLabel', uTL3);
        title('sum of observation potential');
        %
        subplot(1,spSize,4)
        
        imagesc(scoreTemp, [0 1]);
        set(gca,'YDir','normal')
        ylabel('\lambda_1');
        xlabel('\lambda_2');
        % set(gca,'YTick', 1:length(uNL1),'YTickLabel', nL1labels);
        % set(gca,'XTick', 1:length(uTL3),'XTickLabel', uTL3);
        title('score');
        idx = find(score == max(max(score)));
        for i = 1:length(idx)
            idxTemp = idx(i);
            maxTL3 = ceil(idxTemp/size(score,1));
            maxNl1 = mod(idxTemp, size(score,1));
            if maxNl1 == 0
                maxNl1 = size(score,1);
            end
            text(maxTL3-0.2,maxNl1-0.1,'*', 'FontSize', 20)
        end
        
        
    end
end
end

%% makeColor
function cmap = makeCmap()
%%
cSize = 100;
cmap = [];

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
colormap(cmap)
% clear h offset upOffset cmap cSize  cmax

end