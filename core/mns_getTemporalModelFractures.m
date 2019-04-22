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
function mnsResultsSummary = mns_getTemporalModelFractures(mnsResults, metName)
nL1 = zeros(length(mnsResults),1);
nL3 = zeros(length(mnsResults),1);
mns = mnsResults{1};

%make name of RPs
for i = 1:size(mns.results.reactionPairs,1)
    label{i} = [metName{mns.results.reactionPairs(i,1)} ' <-> ' metName{mns.results.reactionPairs(i,2)}];
end

for i = 1:length(mnsResults)
    breakingReaction(:,i) = mnsResults{i}.results.breakingReactions;
    clustLabel(:,i) = mnsResults{i}.results.inferedLabels;
    nL(i) = mnsResults{i}.parameters.nL1;
end



end