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
function [mns, res2] = mns_runTimeFrame(model, dataStruct, initMNS, mnsDataFolder, folder)
% Hardcoded variables
fnProb = 'nCliqueProb.txt';
fObs = 'observedEdges.txt';
fCliques = 'cliques.txt';
fTime = 'temporalEdges.txt';
% alpha = 0.01; %learning rate;
% sigma = 10;
rep = 0;

mnsFolder = regexprep(which('mns_runTimeFrame'),'mns_runTimeFrame\.m','');
mnsBaseFolder = regexprep(which('mns_initialize.m'),'mns_initialize\.m','');
if nargin < 3 || isempty(initMNS)
    initMNS = mns_generateInitMNS();
end
if nargin < 4 || isempty(mnsDataFolder)
    switch computer
        case 'PCWIN32'
            mnsDataFolder = fullfile(mnsFolder, 'data');
        case 'PCWIN64'
            mnsDataFolder = fullfile(mnsFolder, 'data');
        case 'MACI64'
            mnsDataFolder = fullfile(mnsFolder, 'data');
        otherwise
            mnsDataFolder = fullfile(mnsFolder, 'data');
    end
end
% update to mns executable folder
if nargin < 5 || isempty(folder)
    folder = fullfile(mnsBaseFolder, 'opengm');
end
switch computer
    case 'PCWIN32'
        mnsExec = 'mns_timeFrame.exe';
        mnsDataFolder = [mnsDataFolder '\'];
    case 'PCWIN64'
        mnsExec = 'mns_timeFrame_64.exe';
        mnsDataFolder = [mnsDataFolder '\'];
    case 'MACI64'
        mnsExec = 'mns_timeFrame_MAC';
        mnsDataFolder = [mnsDataFolder '/'];
    otherwise
        mnsDataFolder = [mnsDataFolder '/'];
end
solver = initMNS.input.inferenceAlgorithm;

stdType = initMNS.stdType; %options: 'fix', 'one group', 'multiple group'
stdFix = initMNS.stdVal;

normObsProb = initMNS.normObsProb;
neighProbFuncType = initMNS.neighProbFuncType;
biolReplicates = initMNS.biolReplicates;
plotOut = false;

%get Model information from input
iaMat = model.iaMat;
modelMetID = model.metaboliteId;

%get InputData
data = dataStruct.data; %e.g. log_2(fc)
nTp = size(data,2);

%% initalize weights
%nL1 = 1; %weight of the neighbourhood function
%oL2 = 1; %weight of the observation function
%tL3 = 1; %weight of the temporal function
switch initMNS.weights
    case 'uniform'
        nL1 = 10^random('unif', log10(initMNS.nL1min), log10(initMNS.nL1max));
        oL2 = 10^random('unif', log10(initMNS.oL2min), log10(initMNS.oL2max));
        tL3 = 10^random('unif', log10(initMNS.tL3min), log10(initMNS.tL3max));
    case 'exact'
        nL1 = initMNS.nL1min;
        oL2 = initMNS.oL2min;
        tL3 = initMNS.tL3min;
    case 'gaussian'
        nL1 = random('normal', initMNS.nL1min, initMNS.nL1max);
        oL2 = random('normal', initMNS.oL2min, initMNS.oL2max);
        tL3 = random('normal', initMNS.tL3min, initMNS.tL3max);
        while nL1 < 0
            nL1 = random('normal', initMNS.nL1min, initMNS.nL1max);
        end
        while oL2 < 0
            oL2 = random('normal', initMNS.oL2min, initMNS.oL2max);
        end
        while tL3 < 0
            tL3 = random('normal', initMNS.tL3min, initMNS.tL3max);
        end
    otherwise
        
end

%% initalize labeling

switch initMNS.temporalClasses
    case 'global'
        initLabels = zeros(size(data,1)*size(data,2),1);
        switch initMNS.labeling
            case 'kmeans'
                success = false;
                while ~success
                    try
                        initLabels = kmeans(reshape(data,size(data,1)*size(data,2),1), initMNS.noOfclusters, 'replicates', 50)-1;
                        success = true;
                    catch abc
                        disp(abc);
                    end
                end
                
            case 'random'
                initLabels = random('unid',initMNS.noOfclusters, length(initLabels), 1)-1;
            case 'fix'
                    initLabels = initMNS.input.initLabels;
            otherwise
        end
        %calc mean and std of groupings
        [groupMean, groupStd, initLabels] = mns_updateMeanAndStd(reshape(data,size(data,1)*size(data,2),1), initLabels, initMNS.noOfclusters, initMNS.stdType, ...
            biolReplicates, initMNS.stdVal, initMNS.meanType,  initMNS.mean, false);
        groupMean(1:initMNS.noOfclusters,1:size(data,2)) = repmat(groupMean(:),1,size(data,2));
        groupStd(1:initMNS.noOfclusters,1:size(data,2)) = repmat(groupStd(:),1,size(data,2));
        initLabels = reshape(initLabels, size(data,1), size(data,2));
    case 'individual'
        initLabels = zeros(size(data,1),size(data,2));
        groupMean = zeros(initMNS.noOfclusters,size(data,2));
        groupStd = zeros(initMNS.noOfclusters,size(data,2));
        for t = 1:size(data,2)
            switch initMNS.labeling
                case 'kmeans'
                    initLabels(:,t) = kmeans(data(:,t), initMNS.noOfclusters, 'replicates', 50)-1;
                case 'random'
                    initLabels(:,t) = random('unid',initMNS.noOfclusters, length(initLabels), 1)-1;
                case 'fix'
                    initLabels(:,t) = initMNS.input.initLabels(:,t);
                otherwise
            end
            %calc mean and std of groupings
            [groupMean(:,t), groupStd(:,t)] = mns_updateMeanAndStd(data(:,t), initLabels(:,t), initMNS.noOfclusters, initMNS.stdType, ...
                biolReplicates, initMNS.stdVal, initMNS.meanType,  initMNS.mean, false);
        end
        
end
mns.parameters.neighProbFuncType = initMNS.neighProbFuncType;
mns.parameters.initMNS = initMNS;
mns.parameters.noOfLabels = initMNS.noOfclusters;
mns.parameters.nL1 = nL1;
mns.parameters.oL2 = oL2;
mns.parameters.tL3 = tL3;
mns.parameters.biolReplicates = biolReplicates;
mns.parameters.normObsProb = normObsProb;
mns.parameters.groupMean = groupMean;
mns.parameters.groupStd = groupStd;
mns.results.initLabels = initLabels;
mns.results.groupMean = groupMean;
mns.results.groupStd = groupStd;

%% generate cliques from iaMat
if ~isfield(initMNS.input, 'edgeCliqueFileExists') || ~initMNS.input.edgeCliqueFileExists
    mns.nCliques = mns_iaMat2EdgeFile(iaMat, fullfile(mnsDataFolder, fCliques), size(dataStruct.data,2));
else
    mns.nCliques = initMNS.input.nCliques;
end

if isfield(initMNS, 'verbose') && initMNS.verbose
    verbose = ' -v';
else
    verbose = '';
end
% verbose = ' -v';
%calculate potential for neighbouring cliques;
mns.nCliques.Potentials = mns_neighbourProb2nProbFile(fullfile(mnsDataFolder,fnProb), nL1, ...
    mns.nCliques.minClSize, mns.nCliques.maxClSize, neighProbFuncType, true);

%% calculate observation potentials
if ~biolReplicates
    mns.observations = mns_temporalDiffData2ObsEdgeFile(modelMetID, dataStruct.data, dataStruct.annotation, ...
        dataStruct.dataType, groupMean, groupStd, fullfile(mnsDataFolder,fObs), initMNS, oL2, 2, mns.nCliques, ...
        fullfile(mnsDataFolder,fTime), tL3, initLabels);
    mns.observations.obsFacId(mns.observations.obsFacId ~= 0) = mns.observations.obsFacId(mns.observations.obsFacId ~= 0) + max(mns.nCliques.clId);
    mns.observations.biolReplicates = false;
else %with biological replicates
    disp('not yet implemented');
%     mns.observations = mns_diffDataBiolRep2ObsEdgeFile(modelMetID, dataStruct.data, ...
%         dataStuct.annotation, dataStruct.dataType, groupMean, groupStd, [folder fObs], normObsProb, oL2, 2);
%     mns.observations.obsFacId(mns.observations.obsFacId ~= 0) = mns.observations.obsFacId(mns.observations.obsFacId ~= 0) + max(mns.nCliques.clId);
%     mns.observations.biolReplicates = true;
end
% 
% cd(folder)
% dos([folder 'mns_cpp_infer_multipleTimeFrames.exe']);
% cd 'W:\users\Andreas\MATLAB';
% 
if folder == '1'
    mnsExecFullPath = mnsExec;
else
    mnsExecFullPath = fullfile(folder, mnsExec);
end
% tempFolder = regexprep([mnsDataFolder filesep], '\\','\\\\');
tempFolder = mnsDataFolder;
if isfield(initMNS.input, 'inferenceParameter')
    infParameter =  [' -p ' num2str(initMNS.input.inferenceParameter)];
else
    infParameter = '';
end
switch computer
    case 'MACI64'
        unix([mnsExecFullPath ' -i ' solver verbose ' -f ' tempFolder infParameter]);
    otherwise
        dos([mnsExecFullPath ' -i ' solver verbose ' -f ' tempFolder infParameter]);
end
% cd 'W:\users\Andreas\MATLAB';
if exist(fullfile(mnsDataFolder, 'out.txt'),'file')
    mns.results = mns_loadResult2state(mnsDataFolder);
else
    mns.results = [];
    disp(mnsDataFolder);
    return;
end
mns.results.initLabels = initLabels;
[mns.results.breakingReactions, mns.results.breakingReactionTime, mns.results.reactionPairs, mns.results.reactionPairsTime] = ...
    getBreakPointsTemporalModel(model, mns.results.inferedLabels, initMNS.input.reactionPairs,nTp);
mns.results.groupMean = groupMean;
mns.results.groupStd = groupStd;
if plotOut
%     detected = mean(reshape(mns.observations.detected,8,44)',2);
    labels2plot = reshape(mns.results.inferedLabels, length(mns.results.inferedLabels)/size(data,2), size(data,2));
    data2plot = reshape(mns.observations.data, size(data,2),length(mns.results.inferedLabels)/size(data,2))';
    figure, 
    subplot(1,3,1)
    imagesc(data2plot,[-1 1]);
    cmap = redgreencmap(11);
    colormap(gca,cmap);
    title('Data');
    xlabel('Time window');
    ylabel('Metabolites');
%     freezeColors;
    subplot(1,3,2)
    imagesc(labels2plot);
    if initMNS.noOfclusters == 3
        cmap2 = [0.2 0.8 0.8; 0.4 0.4 0.4; 0.8 0.4 0.8];
    else
        cmap2 = hsv(initMNS.noOfclusters);
    end
    colormap(gca,cmap2);
    title(['Labels - \lambda_1 = ' num2str(mns.parameters.nL1) '; \lambda_3 = ' num2str(mns.parameters.tL3)]);
    xlabel('Time window');
    ylabel('Metabolites');
    subplot(1,3,3)
    imagesc(mns.results.breakingReactionTime,[0 1]);
    colormap(gca,[1 1 1; 0.75 0 0]);
    set(gca,'XTickLabels', [1:nTp-1]+0.5)
    xlabel('Time window');
    ylabel('Metabolites');
    title('Temporal Fractures');
    
    
    figure
    for i = 1:size(mns.results.reactionPairs,1)
        if max(mns.results.reactionPairs(i,1)) <= length(model.metaboliteId) && max(mns.results.reactionPairs(i,2)) <= length(model.metaboliteId)
            if isfield(model, 'mat')
                idxRp = find(model.mat(:,mns.results.reactionPairs(i,1)) ~= 0 & model.mat(:,mns.results.reactionPairs(i,2)) ~= 0);
                label{i} = [model.rpId{idxRp} ': ' model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}];
            else
                label{i} = [model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}];
            end
        end
    end
    imagesc(mns.results.breakingReactions, [0 1])
    colormap(gca, [1 1 1; 0.75 0 0])
    title('Fractures')
%     set(gca,'XTick',[1:nTp-1],'XTickLabels', [1:nTp-1]+0.5)
    set(gca,'YTick',1:length(label),'YTickLabel', label)
    xlabel('Time window')
end
% set(gca,'YTick', 1:44)
% set(gca,'YTickLabel', metNames)

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

function [breakingReaction, breakingReactionTime, pairs, pairsTime] ...
    = getBreakPointsTemporalModel(model, labels, pairs, nTp)


if nargin < 3
    pairs = [];
end

if isempty(pairs)
    getPairs = true;
else
    getPairs = false;
end
%% get reaction pairs
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
    %add the pairs for each time point
    pairsTemp = pairs;
    for t = 1:nTp-1
        pairs = [pairs; pairsTemp+t*length(model.metaboliteId)];
    end
end
%% get neighbor breaking reactions for each time windows
breakingReaction = zeros(size(pairs,1),1);
for t = 1:nTp-1
    for i = 1:size(pairs,1)
        if labels(pairs(i,1)) ~= labels(pairs(i,2))
            breakingReaction(i) = 1;
        end
    end
end
breakingReaction = reshape(breakingReaction,length(breakingReaction)/nTp,nTp);
%% make TimePairs

%make pairs
pairsTime = zeros(length(model.metaboliteId), nTp);
for t = 1:nTp
    pairsTime(:,t) = [1:length(model.metaboliteId)]+(t-1)*length(model.metaboliteId);
end

%get temporal breakingReactions
breakingReactionTime = zeros(size(pairsTime,1),nTp-1);
for t = 1:nTp-1
    for i = 1:size(pairsTime,1)
        if labels(pairsTime(i,t)) ~= labels(pairsTime(i,t+1))
            breakingReactionTime(i,t) = 1;
        end
    end
end




end