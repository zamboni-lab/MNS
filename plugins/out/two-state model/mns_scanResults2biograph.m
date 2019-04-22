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
function mns_scanResults2biograph(model, mnsScanResults, sortBy, idx, nNeighborMetabolites, cmax)
if nargin < 3 || isempty(sortBy)
    sortBy = 'rankMax';
end
if nargin < 4 || isempty(idx)
    idx = 1;
end
if nargin < 5 || isempty(nNeighborMetabolites)
    nNeighborMetabolites = 10;
end
if nargin < 6 || isempty(nNeighborMetabolites)
    cmax = 2.5;
end
%%
sPaths = graphallshortestpaths(sparse(model.iaMat), 'Directed' ,0);

%% get sorted results
outTable = mns_scanResult2table(mnsScanResults, model,sortBy);
rpIds = outTable(idx+1,2);
% rpIds = regexp(rpIds, ', ','split');
% rpIds = rpIds{1};
geneSymbols = outTable(idx+1,4);
ranks = outTable(idx+1,1);
metIdx = [];
edgeLabelArr = [];
edgeLabelName = [];
for i = 1:length(rpIds)
    rpIdsTemp = regexp(rpIds{i}, ', ','split');
    rpIdsTemp = rpIdsTemp{1};
    if ischar(rpIdsTemp) == 1
        rpModelIdx(i) = find(strcmp(model.rpId, rpIdsTemp) == 1);
    else
        rpModelIdx(i) = find(strcmp(model.rpId, rpIdsTemp(1)) == 1);
    end
    nMetIdx = find(model.mat(rpModelIdx(i),:) ~= 0);
    nMetIdxArr(i,:) = nMetIdx;
    edgeLabelName = [edgeLabelName; {[model.metaboliteName{nMetIdx(1)} ' -> ' model.metaboliteName{nMetIdx(2)}]}];
    edgeLabelName = [edgeLabelName; {[model.metaboliteName{nMetIdx(2)} ' -> ' model.metaboliteName{nMetIdx(1)}]}];
    edgeLabelArr = [edgeLabelArr; {[num2str(ranks{i}) ': ' geneSymbols{i}]}; {[num2str(ranks{i}) ': ' geneSymbols{i}]}];
    sPathsTemp = min(sPaths(nMetIdx,:));
    temp = sort(sPathsTemp);
    distCo = temp(nNeighborMetabolites);
    if isinf(distCo)
        distCo = max(temp(isinf(temp) == 0));
    end
    metIdx = [metIdx find(sPathsTemp <= distCo)];
end
metIdx = unique(metIdx);

%% reduce iaMat and makeBioGraph
iaMatTemp = model.iaMat(metIdx, metIdx);
metNamesTemp = model.metaboliteName(metIdx);
BGObj2 = biograph(iaMatTemp, metNamesTemp,'LayoutType','equilibrium');

%% make color values
mns = mnsScanResults.mns{1};
fc = mns.observations.data(metIdx);
detected = mns.observations.detected(metIdx);

%% make dualColorColormap
N = 255;
X = [0.5: -1/(N-1):-0.5];
X = abs(X).*2;
cneg = 1-[0 1 0];
cpos = 1-[1 0 0];
R = [1-X(:,1:N/2)*cneg(1) 1-X(:,(N/2 + 1):N)*cpos(1)];
G = [1-X(:,1:N/2)*cneg(2) 1-X(:,(N/2 + 1):N)*cpos(2)];
B = [1-X(:,1:N/2)*cneg(3) 1-X(:,(N/2 + 1):N)*cpos(3)];
cmap = [R' G' B'];
cneg = [0 102 204]/255;
cpos = [204 0 51]/255;

n = 255;

cmap = exp_colormap('custom white', n, cneg, cpos);
%%
% cmap = redgreencmap(255);
% cmax = 2.5;
colorVals = fc;
colorVals(isnan(colorVals)) = min(colorVals) - abs(min(colorVals));
colorVals(isinf(colorVals)) = NaN;
colorVals(isnan(colorVals)) = max(colorVals) + abs(max(colorVals));
maxVal = cmax;
cmpdColorAxisMin = -maxVal;
cmpdColorAxisMax = maxVal;
colors = round((colorVals-cmpdColorAxisMin)/(cmpdColorAxisMax-cmpdColorAxisMin)*253)+1;
colors(colors<1)= 1;
colors(colors>254)= 254;
cmpdColor = cmap(colors,:);
%%
% BGObj2.ArrowSize = 0;
BGObj2.ShowArrows = 'off';
nEdges = length(BGObj2.Edges);
for i = 1:nEdges
    idxTemp = find(strcmp(edgeLabelName, BGObj2.Edges(i).ID) == 1);
    if ~isempty(idxTemp)
        BGObj2.Edges(i).Label = edgeLabelArr{idxTemp};
        BGObj2.Edges(i).LineWidth = 4;
        BGObj2.Edges(i).LineColor = [0.8 0 0];
    else
        BGObj2.Edges(i).Label = '';
    end
end
for i = 1:length(metIdx)
    if detected(i) == 1
        BGObj2.Nodes(i).Color = cmpdColor(i,:);
        BGObj2.Nodes(i).LineColor = [0.2 0.2 0.2];
    else
        BGObj2.Nodes(i).Color = [0.5 0.5 0.5];
        BGObj2.Nodes(i).LineColor = [0.2 0.2 0.2];
    end
%     BGObj2.Nodes(i).Size = [cmpdSize(i) cmpdSize(i)];
end

view(BGObj2)


end