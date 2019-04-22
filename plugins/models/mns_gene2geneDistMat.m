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
function [gene2geneDist,gene2rp] = mns_gene2geneDistMat(model,rpGeneCutOff)
if nargin < 2 || isempty(rpGeneCutOff)
    rpGeneCutOff = inf;
end
gene2geneDist = inf(length(model.GeneName), length(model.GeneName));
met2metDist = graphallshortestpaths(sparse(model.iaMat), 'Directed', false);
gene2rp = zeros(length(model.GeneName), length(model.rpId));
idx2remove = [];
if isfield(model, 'gene2rp')
    idx2remove = find(sum(model.gene2rp,1) >= rpGeneCutOff);
end
for i = 1:length(model.GeneName)
    gene2geneDist(i,i) = 0;
    idxEC = find(model.GeneToEC(i,:) ~= 0);
    idxRP = [];
    idxMetIn = [];
    if ~isempty(idxEC)
        for k = 1:length(idxEC)
            idxRP = [idxRP find(model.ECtoRP(idxEC(k),:) ~= 0)];
            if ~isempty(idxRP)
                gene2rp(i,idxRP) = 1;
                for r = 1:length(idxRP)
                    if isempty(find(idx2remove == idxRP(r)))
                        idxMetIn = [idxMetIn find(model.mat(idxRP(r),:) ~= 0)];
                        
                    end
                end
            end
        end
    end
    if ~isempty(idxMetIn)
        for j = i+1:length(model.GeneName)
            idxEC = find(model.GeneToEC(j,:) ~= 0);
            idxRP = [];
            idxMetOut = [];
            if ~isempty(idxEC)
                for k = 1:length(idxEC)
                    idxRP = [idxRP find(model.ECtoRP(idxEC(k),:) ~= 0)];
                    if ~isempty(idxRP)
                        for r = 1:length(idxRP)
                            if isempty(find(idx2remove == idxRP(r)))
                                idxMetOut = [idxMetOut find(model.mat(idxRP(r),:) ~= 0)];
                            end
                        end
                    end
                    
                end
            end
            if ~isempty(idxMetOut)
                gene2geneDist(i,j) = min(min(met2metDist(idxMetIn, idxMetOut)))+1;
                gene2geneDist(j,i) = gene2geneDist(i,j);
            end
        end
    end
end


end