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
function mnsScanResults = mns_scan2state_multipleParameterizations(model,dataStruct, nameTag, nlRange, varargin)
if isempty(nameTag)
    nameTag = 'temp_mns';
end

if isempty(nlRange)
    nlRange = 0;
end

nParas = length(varargin);
%%
for i = 1:nParas
    mnsScanResults = mns_scan2state(model,dataStruct, [nameTag '_' num2str(i)], nlRange, varargin{i});
    maxBreakingReaction(:,i) = mnsScanResults.maxBreakingReaction;
    sumBreakingReaction(:,i) = mnsScanResults.sumBreakingReaction;
    rankMaxBreaking(:,i) = tiedrank(-mnsScanResults.maxBreakingReaction);
    rankSumBreaking(:,i) = tiedrank(-mnsScanResults.sumBreakingReaction);
    mns{i} = mnsScanResults.mns{1};
end
mnsScanResults.maxBreakingReaction = maxBreakingReaction;
mnsScanResults.sumBreakingReaction = sumBreakingReaction;
mnsScanResults.rankMaxBreaking = rankMaxBreaking;
mnsScanResults.rankSumBreaking = rankSumBreaking;
mnsScanResults.combinedParameterizations = true;
mnsScanResults.initMnsStruct = varargin;
mnsScanResults.mns = mns;
%% calc rankproduct and p-values

nReac = size(rankMaxBreaking,1);
mnsScanResults.rankproductMaxBreaking = prod(rankMaxBreaking,2);
mnsScanResults.rankproductSumBreaking = prod(rankSumBreaking,2);
mnsScanResults.pRankproductMaxBreaking = rankprodbounds(mnsScanResults.rankproductMaxBreaking', nReac, nParas, 'geometric')';
mnsScanResults.pRankproductSumBreaking = rankprodbounds(mnsScanResults.rankproductSumBreaking', nReac, nParas, 'geometric')';
    
end