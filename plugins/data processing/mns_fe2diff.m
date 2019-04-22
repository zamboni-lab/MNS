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
function diffout = mns_fe2diff(pertVec, repVec, ctrl, treat, mat, meanCtrl)

if nargin < 6
    meanCtrl = true;
end

ctrlIdx = [];
treatIdx = [];
rCtrl = regexp(pertVec, ctrl);
rTreat = regexp(pertVec, treat);
for i = 1:length(pertVec)
    if rCtrl{i} == 1
        ctrlIdx = [ctrlIdx; i ];
    elseif rTreat{i} == 1
        treatIdx = [treatIdx; i ];
    end
end
ctrlMat = mat(:,ctrlIdx);
treatMat = mat(:,treatIdx);
ctrlRep = repVec(ctrlIdx);
treatRep = repVec(treatIdx);

uCtrlRep = unique(ctrlRep);
uTreatRep = unique(treatRep);
diffout = zeros(size(mat,1), length(uTreatRep));

for i = 1:length(uTreatRep)
    for j = 1:size(diffout,1)
        if meanCtrl
            if ischar(uTreatRep{i})
                diffout(j,i) = log2(mean(treatMat(j,find(strcmp(treatRep, uTreatRep{i}) == 1)),2)./mean(ctrlMat(j,:),2));
            else
                diffout(j,i) = log2(mean(treatMat(j,treatRep == uTreatRep(i)),2)./mean(ctrlMat(j,:),2));
            end
        else
            diffout(j,i) = log2(mean(treatMat(j,treatRep == uCtrlRep(i)),2)./mean(ctrlMat(j,ctrlRep == uCtrlRep(i)),2));
        end
    end
end


end