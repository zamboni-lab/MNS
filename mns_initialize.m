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
function mns_initialize()
mnsFolder = regexprep(which('mns_initialize.m'),'mns_initialize\.m','');
addpath(fullfile(mnsFolder,'external'));
addpath(fullfile(mnsFolder,'opengm'));
addpath(fullfile(mnsFolder,'core'));
addpath(fullfile(mnsFolder,'models'));
addpath(fullfile(mnsFolder,'examples'));
addpath(fullfile(mnsFolder,'plugins','out'));
addpath(fullfile(mnsFolder,'plugins','data processing'));
addpath(fullfile(mnsFolder,'plugins','models'));
addpath(fullfile(mnsFolder,'plugins','out', 'two-state model'));
addpath(fullfile(mnsFolder,'plugins','out', 'subgraph extraction'));
addpath(fullfile(mnsFolder,'plugins','out', 'continous model'));
addpath(mnsFolder);
end