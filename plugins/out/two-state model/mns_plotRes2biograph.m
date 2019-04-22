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
function mns_plotRes2biograph(mnsResults, model, plotCluster)

if plotCluster
    BGObj2 = biograph(model.iaMat, model.metaboliteName);
    cmap = hsv(mnsResults.parameters.noOfLabels);
    
    for i = 1:length(mnsResults.results.inferedLabels)
        BGObj2.Nodes(i).Color = cmap(mnsResults.results.initLabels(i)+1,:);
    end
    
    view(BGObj2)
end

BGObj = biograph(model.iaMat, model.metaboliteName);
cmap = hsv(mnsResults.parameters.noOfLabels);

for i = 1:length(mnsResults.results.inferedLabels)
    BGObj.Nodes(i).Color = cmap(mnsResults.results.inferedLabels(i)+1,:);
end

view(BGObj)

end