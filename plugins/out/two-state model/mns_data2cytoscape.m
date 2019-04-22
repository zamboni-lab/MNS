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
function mns_data2cytoscape(data, nameBase, nameTag)

outFile = [nameBase ' - input data.xls'];

useTag = true;

if nargin < 3
    nameTag = '';
    useTag = false;
end

if useTag
    outMat(1,:) = {'metId', [nameTag '_log2fc']};
else
    outMat(1,:) = {'metId', 'log2fc'};
end

switch data.dataType
    case 'fiaExp'
        for i = 1:length(data.annotation)
            for j = 1:length(data.annotation(i).id)
                outMat(end+1,1) = data.annotation(i).id(j);
                outMat{end,2} = data.data(i);
            end
        end
    otherwise
        disp('not yet implemented')
end

xlswrite(outFile, outMat);

end