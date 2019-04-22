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
function mns_2stateRes2cytoscape(model, mnsResults, nameBase)

edgeFile = [nameBase ' - mns_breakingReactions.xls'];
nodeFile = [nameBase ' - mns_inferedLabels.xls'];

nMet = length(model.metaboliteId);

metOutMat(:,1) = [{'MetId'}; model.metaboliteId];
pairOutMat = {'rpId'};

mnsRpTemp = mnsResults{1}.results.reactionPairs;
for i = 1:size(mnsResults{1}.results.reactionPairs,1)
    pairOutMat(i+1,1) = {[model.metaboliteId{mnsRpTemp(i,1)} ' - ' model.metaboliteId{mnsRpTemp(i,2)}]}; 
end

for i = 1:length(mnsResults)
    pairOutMat(:,i+1) = [{num2str(i)}; num2cell(mnsResults{i}.results.breakingReactions)];
    metOutMat(:,i+1) = [{num2str(i)}; num2cell(mnsResults{i}.results.inferedLabels)];
end

xlswrite(edgeFile, pairOutMat);
xlswrite(nodeFile, metOutMat);

end