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
function mns_plotTimeCourseLabelOverlay(mns, rawData, rawDataLabels, metNames, plotSingle)

uTime = unique(rawDataLabels);
tempMat = zeros(size(rawData,1),length(uTime));
for i = 1:length(uTime)
    tempMat(:,i) = nanmean(rawData(:,rawDataLabels == uTime(i)),2);
    tempSD(:,i) = nanstd(rawData(:,rawDataLabels == uTime(i)),[],2);
end
if nargin < 5
    plotSingle = false;
end
labels2plot = reshape(mns.results.inferedLabels, ...
    length(mns.results.inferedLabels)/mns.observations.noOfTimePoints, mns.observations.noOfTimePoints);
% data2plot = reshape(mns.observations.data, ...
%     mns.observations.noOfTimePoints,length(mns.results.inferedLabels)/mns.observations.noOfTimePoints)';

detected = mean(reshape(mns.observations.detected,mns.observations.noOfTimePoints,length(mns.results.inferedLabels)/mns.observations.noOfTimePoints),1);
iaId2DataId =  mean(reshape(mns.observations.iaId2dataId,mns.observations.noOfTimePoints,length(mns.results.inferedLabels)/mns.observations.noOfTimePoints),1);
if mns.parameters.noOfLabels == 3
    cmap = [0.2 0.8 0.8; 0.4 0.4 0.4; 0.8 0.4 0.8];
else
    cmap = hsv(mns.parameters.noOfLabels);
end
cmap = [236 0 140; 46 49 146; 0 174 239; 0 166 81;255 242 0]/255;
c = 0;
posF2 = [408.0000  419.0000  140.0000  103.5000];
paperpos = [0,0,10.5,7.7];
posa = [0 0 1 1];
if plotSingle
%     figure
    for i = 1:length(detected)
        if detected(i) == 1
            c = c+1;
%             subplot(5,6,c), hold on
            h = figure;
            hold on
            dataIdx = iaId2DataId(i);
            
            %         colormap(cmap);
            miny = min(rawData(dataIdx,:));
            maxy = max(rawData(dataIdx,:));
            
            if maxy < 0
                maxyBox = maxy*0.8;
                maxyLim = maxy*0.9;
            else
                maxyBox = maxy*1.2;
                maxyLim = maxy*1.1;
            end
            
            if miny < 0
                minyBox = miny*1.2;
                minyLim = miny*1.1;
            else
                minyBox = miny*0.8;
                minyLim = miny*0.9;
            end
            for t = 1:mns.observations.noOfTimePoints
                poly1 = zeros(4,2);
                poly1([1 4],1) = uTime(t);
                poly1([2 3],1) = uTime(t+1);
                poly1([1 2],2) = minyBox;
                poly1([3 4],2) = maxyBox;
                p = patch(poly1(:,1), poly1(:,2), cmap(labels2plot(i,t)+1,:));
                set(p, 'EdgeColor', 'none')
                alpha(p,0.7)
                clear poly1
            end
%             plot(rawDataLabels, rawData(dataIdx,:), '.', 'MarkerSize', 10, 'Color', 'k');
%             plot(uTime, tempMat(dataIdx,:), '-k', 'LineWidth', 3, 'Color', 'k');
            plot(uTime, tempMat(dataIdx,:), '-k', 'LineWidth', 30, 'Color', 'k');
%             xlabel('time','FontSize', 44);
%             ylabel('Ion intensity','FontSize', 44);
            ylim([minyLim maxyLim]);
            xlim([min(uTime) max(uTime)]);
            set(gca, 'FontSize', 6);
            title(metNames(i),'FontSize', 2);
            set(gca, 'XTickLabels', '');
            set(gca, 'YTickLabels', '');
            set(gca, 'Position', posa);
            set(h, 'Position', posF2);
            set(h, 'PaperPosition', paperpos);
            saveas(h, [metNames{i} ' lt=' num2str(mns.parameters.tL3) '.fig'], 'fig');
            saveas(h, [metNames{i} ' lt=' num2str(mns.parameters.tL3) '.png'], 'png');
            close(h);
        end
    end
end

detectedData = labels2plot(detected == 1, :);
figure, hold on 
imagesc(detectedData)
colormap(cmap)
set(gca, 'YTick', 1:size(detectedData,1));
% set(gca, 'YTickLabel', metNames(detected == 1));
xlabel('Timeframes', 'FontSize', 12);
ylim([0.5 size(detectedData,1)+0.5])
xlim([0.5 mns.observations.noOfTimePoints+0.5])
title([{'MNS Results'}; {['\lambda_{n,1} = ' num2str(mns.parameters.nL1) '; \lambda_{o,2} = ' num2str(mns.parameters.oL2) '; \lambda_{t,3} = ' num2str(mns.parameters.tL3) ]}], 'FontSize', 14)

end