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
function mns_plotMultipleScreeningResults(mnsResults, model,wsMax,wnMax,dWs, dWn)
if length(wsMax) == 1
    wtVec = 0:dWs:wsMax;
else
    wtVec = wsMax;
end
if length(wnMax) == 1
    wfVec = 0:dWn:wnMax;
else
    wfVec = wnMax;
end

figure
c = 0;
for i = length(wfVec):-1:1
    for j = 1:length(wtVec)
        c = c+1;
        subplot(length(wfVec), length(wtVec),c)
        hold all
        [prob, score, idx] = ...
            mns_calcProbabilityScanData(mnsResults, model,wtVec(j),wfVec(i),false,false);
        mns_plotTemporalScreeningResults(mnsResults,prob, score,true)
        title(num2str(wtVec(j)))
        ylabel(num2str(wfVec(i)))
        xlabel('')
    end
end

end