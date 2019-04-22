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
function outmodel = rpair2bipartiteModel(model,splitGeneRp)
if nargin < 2
    splitGeneRp = false;
end
%%
if splitGeneRp
     noOfRp = sum(sum(model.gene2rp));
     rpId = cell(1,noOfRp);
     gene2rp = zeros(size(model.gene2rp,1),noOfRp);
     mat = zeros(noOfRp, length(model.metaboliteId));
     rpC = 0;
     for i = 1:length(model.rpId)
         idx = find(model.gene2rp(:,i) ~= 0);
         if ~isempty(idx)
             for j = 1:length(idx)
                 rpC = rpC+1;
                 gene2rp(idx(j),rpC) = 1;
                 mat(rpC,:) = model.mat(i,:);
                 rpId{rpC} = sprintf('%s_%03d', model.rpId{i}, j);
             end
         else
             rpC = rpC+1;
             gene2rp(:,rpC) = model.gene2rp(:,i);
             mat(rpC,:) = model.mat(i,:);
             rpId(rpC) = model.rpId(i);
         end
     end
     model.rpId = rpId;
     model.gene2rp = gene2rp;
     model.mat = mat;
end
%%
iaMat = zeros(length(model.rpId)+length(model.metaboliteId));
iaMatId = [model.rpId';model.metaboliteId];
noRp = length(model.rpId);
isRp = [ones(noRp,1);zeros(length(model.metaboliteId),1)];
isMet = [zeros(noRp,1);ones(length(model.metaboliteId),1)];
for i = 1:noRp
    idx = find(model.mat(i,:) ~= 0);
    iaMat(i,idx+noRp) = 1;
    iaMat(idx+noRp,i) = 1;
end
outmodel = model;
outmodel.iaMat = sparse(iaMat);
outmodel.iaMatId = iaMatId;
outmodel.isRp = isRp;
outmodel.isMet = isMet;
end