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
function mns_plot2stateOutOnKegg(labels, gr, cmap, saveFile,addData)
save = true;
if nargin < 4
    save = false;
end
if nargin < 5
    addData = [];
end
    
gr = gr +1;
org = 'hsa';
dotSize = 60;
%start to make output file
output = {['Organism:' org ]}; % Organism ->
output = [output, {['compound,#E0E0E0,' num2str(dotSize) '']}]; % Make all compounds grey
output = [output, {'reaction,#E0E0E0,10'}]; % Make all reactions grey

for iID = 1:length(labels)
    if regexp(labels{iID},'C\d\d\d\d\d')
            output = [output, {sprintf('%s,#%02X%02X%02X,%0.0f',labels{iID},round(cmap(gr(iID),1)*255),...
                round(cmap(gr(iID),2)*255),round(cmap(gr(iID),3)*255),dotSize)}];
    end
end
output = [output addData'];
ppObj = PathwayProjectorService;
[unused,addr] = mapping(ppObj,output);
% if mode.show
    web(addr, '-browser');
% end
savePtwPro(addr, 'full', saveFile)

end