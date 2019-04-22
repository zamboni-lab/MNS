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
function out = mns_iaMat2EdgeFile(iaMat, filename, nTp)

if nargin < 3
    nTp = 1;
end

% identify max cliques of a graph
uIA = triu(iaMat);
maxCl = maximalCliques(uIA);
out.clId = 0:length(maxCl)-1;
% identify max clique size
maxClSize = 0;
minClSize = 1000;
for i = 1:length(maxCl)
    maxClSize = max(length(maxCl{i}), maxClSize);
    minClSize = min(length(maxCl{i}), minClSize);
end
out.minClSize = minClSize;
out.maxClSize = maxClSize;
out.cliqueMat = nan(length(maxCl),maxClSize);

% print cliques to file
fid = fopen(filename,'w');
out.varId(:,1) = 1:size(iaMat,1);
out.varId = out.varId-1;

out.varId2clId = nan(length(out.varId),100);
out.nCliquePerVar = zeros(length(out.varId),nTp);
out.nNeighboorMetabolites = repmat(sum(iaMat,2),1,nTp);
offSet = max(out.varId)+1;
for t = 1:nTp
    out.varId(:,t) = out.varId(:,1)+(t-1)*offSet;
    for i = 0:1:size(maxCl,1)
        
        if i == 0
            if t == 1
                fprintf(fid,'%d\t%d\r\n', size(maxCl,1)*nTp, maxClSize);
            end
        else
            if t > 1
                maxCl{i} = maxCl{i}+offSet;
            end
            for j = 0:length(maxCl{i})
                if j == 0
                    fprintf(fid,'%d', length(maxCl{i}));
                    out.cliqueSize(i,t) = length(maxCl{i});
                else
                    fprintf(fid,'\t%d', maxCl{i}(j)-1);
                    out.clique{i,t}(j) = maxCl{i}(j)-1;
                    out.cliqueMat(i,j,t) = maxCl{i}(j)-1;
                    idxTemp = find(out.varId(:,t) == maxCl{i}(j)-1);
                    if ~isempty(idxTemp)
                        out.nCliquePerVar(idxTemp,t) = out.nCliquePerVar(idxTemp) + 1;
                        out.varId2clId(idxTemp,out.nCliquePerVar(idxTemp,t),t) = i;
                    end
                end
            end
            if ~(i == size(maxCl,1) && t == nTp)
                fprintf(fid,'\r\n');
            end
        end
    end
end
fclose(fid);
out.varId2clId(:,max(out.nCliquePerVar)+1:end) = [];

end