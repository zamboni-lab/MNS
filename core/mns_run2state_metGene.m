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
function [mns, res2, t] = mns_run2state_metGene(model, dataStruct, initMNS, mnsDataFolder, folder)
tic
%get mns base folder
mnsBaseFolder = regexprep(which('mns_initialize.m'),'mns_initialize\.m','');

fnProb = 'nCliqueProb.txt';
fObs = 'observedEdges.txt';
fCliques = 'cliques.txt';

if nargin < 3 || isempty(initMNS)
    initMNS = mns_generateInitMNS();
end
if nargin < 4 || isempty(mnsDataFolder)
    mnsDataFolder = fullfile(mnsBaseFolder,'data');
end
% update to mns executable folder
if nargin < 5 || isempty(folder)
    folder = fullfile(mnsBaseFolder, 'opengm');
    switch computer
        case 'PCWIN32'
            mnsExec = 'mns_2state.exe';
        case 'PCWIN64'
            mnsExec = 'mns_2state_64.exe';
        case 'MACI64'
            mnsExec = 'mns_2state_MAC';
        otherwise
    end
end
solver = initMNS.input.inferenceAlgorithm;
switch computer
    case 'PCWIN32'
        mnsExec = 'mns_2state.exe';
        mnsDataFolder = [mnsDataFolder '\'];
    case 'PCWIN64'
        mnsExec = 'mns_2state_64.exe';
        mnsDataFolder = [mnsDataFolder '\'];
    case 'MACI64'
        mnsExec = 'mns_2state_MAC';
        mnsDataFolder = [mnsDataFolder '/'];
    otherwise
        mnsDataFolder = [mnsDataFolder '/'];
end
normObsProb = initMNS.normObsProb;
neighProbFuncType = initMNS.neighProbFuncType;
biolReplicates = initMNS.biolReplicates;

%% get Model information from input
iaMat = model.iaMat;
modelIaMatID = model.iaMatId;
modelIsRp = model.isRp;
modelIsMet = model.isMet;
modelMetID = model.metaboliteId;
modelRpID = model.rpId;
modelGeneSymbol = model.GeneSymbol;

%% get InputData
data = dataStruct.data; %e.g. log_2(fc)
geneData = dataStruct.geneData; %e.g. log_2(fc)

%% initalize weights
%nL1 = 1; %weight of the neighbourhood function
%oL2 = 1; %weight of the observation function
switch initMNS.weights
    case 'uniform'
        nL1 = 10^random('unif', log10(initMNS.nL1min), log10(initMNS.nL1max));
        oL2 = 10^random('unif', log10(initMNS.oL2min), log10(initMNS.oL2max));
    case 'exact'
        nL1 = initMNS.nL1min;
        oL2 = initMNS.oL2min;
    case 'gaussian'
        nL1 = random('normal', initMNS.nL1min, initMNS.nL1max);
        oL2 = random('normal', initMNS.oL2min, initMNS.oL2max);
        while nL1 < 0
            nL1 = random('normal', initMNS.nL1min, initMNS.nL1max);
        end
        while oL2 < 0
            oL2 = random('normal', initMNS.oL2min, initMNS.oL2max);
        end
    otherwise
        
end

%% initalize labeling
initLabels = zeros(size(data,1),1);
switch initMNS.labeling
    case 'kmeans'
        success = false;
        while ~success
            try
                %% for metabolites
                initLabels = kmeans(data, initMNS.noOfclusters, 'replicates', 50)-1;
                clMean = zeros(initMNS.noOfclusters,1);
                uLabels = 0:initMNS.noOfclusters-1;
                for i = 1:initMNS.noOfclusters
                    clMean(i) = mean(data(initLabels == i-1));
                end
                [~,idxS] = sort(clMean);
                initLabelsTemp = initLabels;
                for i = 1:initMNS.noOfclusters
                    initLabels(initLabelsTemp == idxS(i)-1) = i-1;
                end
                %% for genes
                initLabelsGene = kmeans(geneData, initMNS.noOfclusters, 'replicates', 50)-1;
                clMeanGene = zeros(initMNS.noOfclusters,1);
                uLabelsGene = 0:initMNS.noOfclusters-1;
                for i = 1:initMNS.noOfclusters
                    clMeanGene(i) = mean(geneData(initLabels == i-1));
                end
                [~,idxS] = sort(clMeanGene);
                initLabelsTemp = initLabelsGene;
                for i = 1:initMNS.noOfclusters
                    initLabelsGene(initLabelsTemp == idxS(i)-1) = i-1;
                end
                success = true;
            catch abc
                disp(abc);
            end
        end
    case 'random'
        initLabels = random('unid',initMNS.noOfclusters, size(data,1), 1)-1;
        initLabelsGene = random('unid',initMNS.noOfclusters, size(geneData,1), 1)-1;
    case 'fix'
        initLabels = initMNS.input.initLabels;
        initLabelsGene = initMNS.input.initLabelsGene;
    otherwise
        
end
mns.parameters.neighProbFuncType = neighProbFuncType;
mns.parameters.initMNS = initMNS;
mns.parameters.noOfLabels = initMNS.noOfclusters;
mns.parameters.nL1 = nL1;
mns.parameters.oL2 = oL2;
mns.parameters.biolReplicates = biolReplicates;
mns.parameters.normObsProb = normObsProb;

%% calc mean and std of groupings
% for metabolites
[groupMean, groupStd] = mns_updateMeanAndStd(data, initLabels, initMNS.noOfclusters,...
    initMNS.stdType, biolReplicates, initMNS.stdVal, initMNS.meanType, initMNS.mean, false);
% for genes
[groupGeneMean, groupGeneStd] = mns_updateMeanAndStd(geneData, initLabelsGene, initMNS.noOfclusters, ...
    initMNS.stdType, biolReplicates, initMNS.subGraphExtraction.geneStdVal, ...
    initMNS.meanType, initMNS.subGraphExtraction.geneMean, false);

%% save initial label files into
% mns_initLabels2labelFile(initLabels, [folder fLabels]);

%% generate nighbouring cliques from iaMat if it does not exist
if ~isfield(initMNS.input, 'edgeCliqueFileExists') || ~initMNS.input.edgeCliqueFileExists
    mns.nCliques = mns_iaMat2EdgeFile(iaMat, fullfile(mnsDataFolder, fCliques));
else
    mns.nCliques = initMNS.input.nCliques;
end

if isfield(initMNS, 'verbose') && initMNS.verbose
    verbose = ' -v';
else
    verbose = '';
end
    

%% calculate potential for neighbouring cliques;
mns.nCliques.Potentials = mns_neighbourProb2nProbFile(fullfile(mnsDataFolder, fnProb), nL1, ...
    mns.nCliques.minClSize, mns.nCliques.maxClSize, neighProbFuncType, true);

%% calculate observation potentials
if ~biolReplicates
    %get metabolite observations but don't print to file
    mns.observations = mns_diffData2ObsEdgeFile(modelMetID, dataStruct, groupMean, groupStd, ...
        fullfile(mnsDataFolder, fObs), initMNS, nL1, oL2, false, mns.nCliques, initLabels);
    mns.observations = mns_geneMetDiffData2ObsEdgeFile(dataStruct, model, groupGeneMean, ...
    groupGeneStd, fullfile(mnsDataFolder, fObs), initMNS, nL1, oL2, mns.nCliques, initLabelsGene, mns.observations, true);
    mns.observations.obsFacId(mns.observations.obsFacId ~= 0) = mns.observations.obsFacId(mns.observations.obsFacId ~= 0) + max(mns.nCliques.clId);
    mns.observations.biolReplicates = false;
else %with biological observations
    mns.observations = mns_diffDataBiolRep2ObsEdgeFile(modelMetID, dataStruct.data, ...
        dataStruct.annotation, dataStruct.dataType, groupMean, groupStd, fullfile(mnsDataFolder,fObs),...
        normObsProb, oL2, 2, mns.nCliques, initMNS.biolReplicatesAvg);
    mns.observations.obsFacId(mns.observations.obsFacId ~= 0) = mns.observations.obsFacId(mns.observations.obsFacId ~= 0) + max(mns.nCliques.clId);
    mns.observations.biolReplicates = true;
end
%%
if folder == '1'
    mnsExecFullPath = mnsExec;
else
    mnsExecFullPath = fullfile(folder, mnsExec);
end
%%
if isfield(initMNS.input, 'inferenceParameter')
    infParameter = [' -p '  num2str(initMNS.input.inferenceParameter)];
else
    infParameter = '';
end
switch computer
    case 'MACI64'
        unix([mnsExecFullPath ' -i ' solver verbose ' -f ' mnsDataFolder infParameter]);
    otherwise
        dos([mnsExecFullPath ' -i ' solver verbose ' -f ' regexprep(mnsDataFolder, '\\','\\\\') infParameter]);
end
if exist(fullfile(mnsDataFolder, 'out.txt'),'file')
    mns.results = mns_loadResult2state(mnsDataFolder);
else
    mns.results = [];
    disp(mnsDataFolder);
    return;
end
mns.results = mns_loadResult2state(mnsDataFolder);
mns.results.initLabels = initLabels;
[mns.results.breakingReactions, mns.results.reactionPairs] = ...
    getBreakPoints(model, mns.results.inferedLabels, initMNS.input.reactionPairs);
mns.results.groupMean = groupMean;
mns.results.groupStd = groupStd;

res2 = 0;
t = toc;
end

function printMNSsolution(res)
    figure, hold on
    subplot(2,2,1)
    plot(0:length(res.solutions)-1, res.solutions, 'LineWidth', 2)
    ylabel('Solution', 'FontSize', 12);
    xlabel('Iteration', 'FontSize', 12);
    title('Solution', 'FontSize', 14);
    set(gca,'FontSize', 12);
   
    subplot(2,2,2)
    imagesc(res.labels, [min(min(res.labels)) max(max(res.labels))]);
    cmap = hsv(max(max(res.labels))+1);
    cmap(2:end+1,:) = cmap;
    cmap(1,:) = [0.7 0.7 0.7];
    colormap(cmap)
    ylabel('Metabolite', 'FontSize', 12);
    xlabel('Iteration', 'FontSize', 12);
    title('Solution', 'FontSize', 14);
    set(gca,'FontSize', 12);
    
    subplot(2,2,3)
    plot(0:length(res.nL1list)-1, res.nL1list, 'LineWidth', 2)
    ylabel('\lambda_1', 'FontSize', 12);
    xlabel('Iteration', 'FontSize', 12);
    title('Neighbour Weight', 'FontSize', 14);
    set(gca,'FontSize', 12);
    
    subplot(2,2,4)
    plot(0:length(res.oL2list)-1, res.oL2list, 'LineWidth', 2)
    ylabel('\lambda_2', 'FontSize', 12);
    xlabel('Iteration', 'FontSize', 12);
    title('Observation Weight', 'FontSize', 14);
    set(gca,'FontSize', 12);

end

% checks if labels between between neighbouring metabolites are the same or
% not
% Input:
% model -> metabolic model
% labels -> cluster/module labels
% pairs -> ids of reaction pairs (if empty it gets the rps from the model)
function [breakingReaction, pairs] = getBreakPoints(model, labels, pairs)


if nargin < 3
    pairs = [];
end

if isempty(pairs)
    getPairs = true;
else
    getPairs = false;
end

if getPairs
    for i = 1:size(model.iaMat,1)
        idx = find(model.iaMat(i,:) ~= 0);
        for k = 1:length(idx)
            
            if ~isempty(pairs)
                if isempty(find(ismember(pairs(:,1), idx(k)) == 1 & ismember(pairs(:,2), i) == 1))
                    pairs = [pairs; i, idx(k)];
                end
            else
                pairs = [pairs; i, idx(k)];
            end
        end
    end
end

breakingReaction = zeros(size(pairs,1),1);

for i = 1:size(pairs,1)
    if labels(pairs(i,1)) ~= labels(pairs(i,2))
        breakingReaction(i) = 1;
    end
end

end