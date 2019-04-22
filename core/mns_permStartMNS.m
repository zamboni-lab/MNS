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
function [breakingReactions, noOfBreakingPermutations] = mns_permStartMNS(model, initMNS, dataStruct, mnsFolder, mnsSubDataFolder, nl,mnsExecFolderP)
    if nargin > 5 && ~isempty(nl)
        initMNS.nL1max = nl;
        initMNS.nL1min = nl;
    end
    
    mkdir(mnsSubDataFolder);
%     disp(fullfile(mnsFolder, 'nCliqueNoPerm.txt'))
    copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'), fullfile(mnsSubDataFolder, 'nCliqueNoPerm.txt'));
    if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
        copyfile(initMNS.input.edgeCliqueFile, fullfile(mnsSubDataFolder ,'cliques.txt'));
    end
    [mnsTEMP, ~, ~] = mns_run2state(model, dataStruct, initMNS, mnsSubDataFolder,mnsExecFolderP);
    breakingReactions = mnsTEMP.results.breakingReactions;
    noOfBreakingPermutations = sum(mnsTEMP.results.breakingReactions);
    try
        rmdir(mnsSubDataFolder,'s');
    catch err
        disp(err)
    end
end




