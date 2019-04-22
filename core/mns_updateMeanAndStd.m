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
function [groupMean, groupStd, groups] = mns_updateMeanAndStd(data, groups, noOfclusters, stdType, biolReplicates, stdFix, meanType, meanVal, plotGrouping)


groupMeanTemp = zeros(1,noOfclusters);
groupStd = zeros(1,noOfclusters);
if ~biolReplicates
    for i = 0:noOfclusters-1
        groupMeanTemp(i+1) =  mean(data(groups == i,:));
        idxTemp(i+1) = i;
    end
    [groupMean, idxS] = sort(groupMeanTemp);
    groupTemp = groups;
    idx = idxTemp(idxS);
    for i = 0:noOfclusters-1
        groups(groupTemp == idx(i+1)) = i;
    end
%     disp(stdType)
    switch stdType
        case 'one group'
            tempVec = zeros(size(data,1),1);
            for i = 0:noOfclusters-1;
                tempVec(groups == i) = groupMean(i+1);
            end
            groupStd(1:end) = sqrt(sum((data-tempVec).^2)/size(data,1));
        case 'multiple groups'
            for i = 0:noOfclusters-1
                groupStd(i+1) =  std(data(groups == i,:));
            end
        case 'fix'
            groupStd(1:end) = stdFix;
        case 'one group - factor'
            tempVec = zeros(size(data,1),1);
            for i = 0:noOfclusters-1;
                tempVec(groups == i) = groupMean(i+1);
            end
            groupStd(1:end) = sqrt(sum((data-tempVec).^2)/size(data,1));
            groupStd  = groupStd/stdFix;
        case 'multiple groups - factor'
            for i = 0:noOfclusters-1
                groupStd(i+1) =  std(data(groups == i,:));
            end
            groupStd = groupStd/stdFix;
        case 'one group - all data'
            groupStd(1:end) = std(data);
%             disp(groupStd)
        case 'one group - all data - factor'
            groupStd(1:end) = std(data)/stdFix;
%             disp(groupStd)
        otherwise
    end
else
    for i = 0:noOfclusters-1
        groupMean(i+1) =  mean(mean(data(groups == i,:)));
    end
    switch stdType
        case 'one group'
            tempVec = zeros(size(data,1),1);
            for i = 1:noOfclusters;
                tempVec(groups == i) = groupMean(i+1);
            end
            tempVec = tempVec';
            for i = 1:size(data,2)
                tempVec(:,i) = tempVec(:,1);
            end
            groupStd(1:end) = sqrt(sum((data-tempVec).^2)/size(data,1));
        case 'multiple groups'
            for i = 0:nOfclusters-1
                groupStd(i+1) =  std(std(data(groups == i,:)));
            end
        case 'fix'
            groupStd(1:end) = stdFix;
        case 'one group - factor'
            tempVec = zeros(size(data,1),1);
            for i = 1:noOfclusters;
                tempVec(groups == i) = groupMean(i+1);
            end
            tempVec = tempVec';
            for i = 1:size(data,2)
                tempVec(:,i) = tempVec(:,1);
            end
            groupStd(1:end) = sqrt(sum((data-tempVec).^2)/size(data,1))/stdFix;
        case 'multiple groups - factor'
            for i = 0:nOfclusters-1
                groupStd(i+1) =  std(std(data(groups == i,:)));
            end
            groupStd = groupStd/stdFix;
        case 'one group - data - factor'
            groupStd(1:end) = std(data)/stdFix;
        otherwise
    end
end
groupStd(groupStd == 0) = eps;
if noOfclusters == 3
%     cmap = redgreencmap(noOfclusters);
    cmap = [0.2 0.8 0.8; 0.4 0.4 0.4; 0.8 0.4 0.8];
else
    cmap = hsv(noOfclusters);
end

switch meanType
    case 'initLabels'
        %do nothing
    case 'fix'
        if length(meanVal) == length(groupMean)
            groupMean = meanVal;
        else
            disp('ERROR: length of fix mean values ~= no of clusters. groupMeans are calculated from initial labels')
        end
    case 'linear - std'
        stepSize = ((max(data)-mean(groupStd))-(min(data)+mean(groupStd)))/(noOfclusters-1);
        groupMean = (min(data)+mean(groupStd)):stepSize:(max(data)-mean(groupStd));
    case 'linear - quantile'
        if isempty(meanVal)
            meanVal = 0.001;
        end
        stepSize = (quantile(data,1-meanVal)-quantile(data,meanVal))/(noOfclusters-1);
        groupMean = quantile(data,meanVal):stepSize:quantile(data,1-meanVal);
end
        
% groupMean = [-1.0 0 1 2 3];
% groupMean = [-1 1 3];
% groupMean = [-2.5:5.5/4:3];
if plotGrouping
    figure, hold on
    xVec = -10:0.01:10;
    for i = 0:noOfclusters-1
        [freq,x] = hist(data(groups == i),10);
        prob = exp( -(xVec-groupMean(i+1)).^2 / (2*groupStd(i+1)^2));
        plot(xVec, prob, '-','LineWidth', 2, 'Color', cmap(i+1,:));
        bar(x, freq/max(freq), 1, 'FaceColor', 'none', 'EdgeColor', cmap(i+1,:));
        legLabels{(i)*2+1} = ['Group ' num2str(i)];
        legLabels{(i)*2+2} = ['Group ' num2str(i) ': Data'];
    end
    xlim([-5 5]);
    legend(legLabels);
end

end
