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
function mnsSubGraphResults = mns_extractSubGraphs_MetGeneIntegration(model,dataStruct, nameTag, nlRange, initMNS, plotResults)
t1 = tic;
%% check if data is logarithmic
if max(dataStruct.data) <= 1 && min(dataStruct.data) >= 0
    dataStruct.data = -log10(dataStruct.data);
end
if max(dataStruct.geneData) <= 1 && min(dataStruct.geneData) >= 0
    dataStruct.geneData = -log10(dataStruct.geneData);
end
%% check input arguments
if nargin < 3 || isempty(nameTag)
    nameTag = 'temp';
end
if nargin < 4 || isempty(nlRange)
    nlRange = 0;
end
if nargin < 5 || isempty(initMNS)
    temp.minData = min(dataStruct.data);
    temp.maxData = max(dataStruct.data);
    temp.std0 = 3/3; 
    temp.std1 = (max(dataStruct.data)-1)/3;
    initMNS = mns_generateInitMNS('modeSubGraphExtraction', temp);
    initMNS.nL1steps = 20;
    initMNS.determinePvalue = 2;
    initMNS.parallel.useCluster = 1;
%     initMNS.permutations = 100;
end
if nargin < 6 || isempty(plotResults)
    plotResults = false;
end
if isfield(dataStruct, 'dataReplicates')
    initMNS.biolReplicates = dataStruct.dataReplicates;
end
%%
% remove all white spaces from the nametag
nameTag = regexprep(nameTag, '\s', '\_');
% get mns folder and create data Folder
mnsFolder = regexprep(which('mns_scan2state.m'),'mns_scan2state\.m','');
mnsBaseFolder = regexprep(which('mns_initialize.m'),'mns_initialize\.m','');
switch computer
    case 'PCWIN32'
        mnsExec = 'mns_2state.exe';
    case 'PCWIN64'
        mnsExec = 'mns_2state_64.exe';
    case 'MACI64'
        mnsExec = 'mns_2state_MAC';
    otherwise
end
if ~isempty(initMNS.input.dataFolder) && isfield(initMNS.input, 'dataFolder')
    mnsDataFolderStem = initMNS.input.dataFolder;
    clFolder = mnsDataFolderStem;
    if ~isempty(initMNS.input.mnsExecFolder) && isfield(initMNS.input, 'mnsExecFolder')
        mnsExecFolder = initMNS.input.mnsExecFolder;
    else
        mnsExecFolder = '';
    end
else
    clFolder = mnsFolder;
    mnsDataFolderStem = fullfile(mnsBaseFolder, 'data');
    mnsExecFolder = '';
end
dataFolders = dir(mnsDataFolderStem);
dataFolders = [{dataFolders(:).name};]';
dataFoldersStem = regexprep(dataFolders, '_\d+$','');
idxSameStem = find(strcmp(dataFoldersStem, nameTag) == 1);
if ~isempty(idxSameStem)
    maxIteration = 0;
    for i = 1:length(idxSameStem)
        a = regexp(dataFolders(idxSameStem(i)), '_(\d+)$','tokens');
        if ~isempty(a{1});
            maxIteration = max(maxIteration, str2num(a{1}{1}{1}));
        end
    end
    nameTag = [nameTag '_' num2str(maxIteration + 1)];
end
switch computer
    case 'PCWIN32'
        mnsDataFolder = [mnsDataFolderStem '\' nameTag '\'];
    case 'PCWIN64'
        mnsDataFolder = [mnsDataFolderStem '\' nameTag '\'];
    case 'MACI64'
        mnsDataFolder = [mnsDataFolderStem '/' nameTag '/'];
    otherwise
        mnsDataFolder = [mnsDataFolderStem '/' nameTag '/'];
end

mkdir(mnsDataFolder)

%% if nlRange = 0, do coarse grained search
if nlRange == 0
    if(initMNS.verboseScan == 2)
        disp('Start coarse grained search...');
    end
    s = 1;
    initMNS.stdVal = [initMNS.subGraphExtraction.std0Steps(s) ...
        (initMNS.subGraphExtraction.maxData - initMNS.subGraphExtraction.std1Steps(s))/3 ...
        (initMNS.subGraphExtraction.maxData - initMNS.subGraphExtraction.std1Steps(s))/3];
    initMNS.subGraphExtraction.geneStdVal = [initMNS.subGraphExtraction.geneStd0Steps(s) ...
        (initMNS.subGraphExtraction.geneMaxData - initMNS.subGraphExtraction.geneStd1Steps(s))/3 ...
        (initMNS.subGraphExtraction.geneMaxData - initMNS.subGraphExtraction.geneStd1Steps(s))/3];
    nlTemp = 1;
    initMNStemp = initMNS;
    i = 0;
    while(nlTemp <= 512)
        i = i+1;
        mnsSubDataFolder = fullfile(mnsDataFolder, ['coarse_nl_' num2str(nlTemp)]);
        mkdir(mnsSubDataFolder);
        copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(mnsSubDataFolder, 'nCliqueNoPerm.txt'));
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, fullfile(mnsSubDataFolder, 'cliques.txt'));
        end
        if(initMNS.verboseScan == 2)
            disp(nlTemp);
        end
        initMNStemp.nL1max = nlTemp;
        initMNStemp.nL1min= nlTemp;
        [mns] = mns_run2state_metGene(model, dataStruct, initMNStemp, mnsSubDataFolder, mnsExecFolder);
        if isempty(mns.results)
            mnsSubGraphResults = [];
            return;
        end
        if i == 1
            initMNStemp.input.initLabels = mns.results.initLabels;
            initMNStemp.labeling = 'fix';
            initMNStemp.input.reactionPairs = mns.results.reactionPairs;
            initMNStemp.input.edgeCliqueFile = fullfile(mnsSubDataFolder, 'cliques.txt');
            initMNStemp.input.edgeCliqueFileExists = true;
            initMNStemp.input.nCliques = mns.nCliques;
            initMNS.input.edgeCliqueFile = fullfile(mnsSubDataFolder, 'cliques.txt');
            initMNS.input.edgeCliqueFileExists = true;
            initMNS.input.nCliques = mns.nCliques;
        end
        nlTemp = nlTemp*2;
        if nlTemp >= 512
            nlRange = 0:nlTemp/initMNS.nL1steps:nlTemp;
            if(initMNS.verboseScan == 2)
                disp('Finished coarsed grain search');
                disp(['nlRange: 0:' num2str(nlTemp/initMNS.nL1steps) ':' num2str(nlTemp)]);
            end
            break;
        end
        
        if sum(mns.results.breakingReactions,1) == 0
            nlRange = 0:nlTemp/initMNS.nL1steps:nlTemp;
            if(initMNS.verboseScan == 2)
                disp('Finished coarsed grain search');
                disp(['nlRange: 0:' num2str(nlTemp/initMNS.nL1steps) ':' num2str(nlTemp)]);
            end
            break;
        end
        
    end
    clear mnsScanResults initMNStemp i mns;
    if(initMNS.verboseScan == 2)
        disp('Coarse grained search finished')
    end
end

%% scanning mode
if(initMNS.verboseScan == 2)
    disp('Start scanning...')
end
subGraphList = [];
for s = 1:length(initMNS.subGraphExtraction.std0Steps) % scan through std values
    initMNS.stdVal = [initMNS.subGraphExtraction.std0Steps(s) ...
        (initMNS.subGraphExtraction.maxData - initMNS.subGraphExtraction.std1Steps(s))/3 ...
        (initMNS.subGraphExtraction.maxData - initMNS.subGraphExtraction.std1Steps(s))/3];
    initMNS.subGraphExtraction.geneStdVal = [initMNS.subGraphExtraction.geneStd0Steps(s) ...
        (initMNS.subGraphExtraction.geneMaxData - initMNS.subGraphExtraction.geneStd1Steps(s))/3 ...
        (initMNS.subGraphExtraction.geneMaxData - initMNS.subGraphExtraction.geneStd1Steps(s))/3];
    for i = 1:length(nlRange)
        mnsSubDataFolder = fullfile(mnsDataFolder, ['std_' num2str(s) '_nl_' num2str(nlRange(i))]);
        mkdir(mnsSubDataFolder);
        copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(mnsSubDataFolder, 'nCliqueNoPerm.txt'));
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, fullfile(mnsSubDataFolder, 'cliques.txt'));
        end
        if(initMNS.verboseScan == 2)
            disp('');
            disp(['Starting step ' num2str((s-1)*length(nlRange)+i) ...
                ', std0 = ' num2str(initMNS.stdVal(1)) ', std1 = ' num2str(initMNS.stdVal(2)) ', nl = ' num2str(nlRange(i))]);
            disp('');
        end
        initMNS.nL1max = nlRange(i);
        initMNS.nL1min= nlRange(i);
        [mns] = mns_run2state_metGene(model, dataStruct, initMNS, mnsSubDataFolder, mnsExecFolder);
        if isempty(mns.results)
            mnsSubGraphResults = [];
            return;
        end
        mnsSubGraphResults.mns(i) = {mns};
        mnsSubGraphResults.breakingReaction(:,i) = mns.results.breakingReactions;
        mnsSubGraphResults.clustLabel(:,i) = mns.results.inferedLabels;
        [subGraphs, listTemp] = getSubgraphs(mns.results.breakingReactions, mns.results.reactionPairs, model.iaMat, mns.observations.iaId2dataId, mns.results.inferedLabels, model.isRp);
        if ~isempty(listTemp)
            subGraphList = [subGraphList; num2cell([ones(size(listTemp,1),1).*initMNS.stdVal(1) ...
                ones(size(listTemp,1),1).*initMNS.stdVal(2)]) num2cell(ones(size(listTemp,1),1).*nlRange(i)) listTemp];
        end
        mnsSubGraphResults.subGraphs(i) = {subGraphs};
        if i == 1
            initMNS.input.edgeCliqueFile = fullfile(mnsSubDataFolder, 'cliques.txt');
            initMNS.input.edgeCliqueFileExists = true;
            initMNS.input.nCliques = mns.nCliques;
            initMNS.input.initLabels = mnsSubGraphResults.mns{1}.results.initLabels;
            initMNS.labeling = 'fix';
            initMNS.input.reactionPairs = mnsSubGraphResults.mns{1}.results.reactionPairs;
            mnsSubGraphResults.maxBreakingReaction = zeros(size(mnsSubGraphResults.breakingReaction,1),1);
            mnsSubGraphResults.sumBreakingReaction = zeros(size(mnsSubGraphResults.breakingReaction,1),1);
        end
        mnsSubGraphResults.maxBreakingReaction(mnsSubGraphResults.breakingReaction(:,i) == 1) = i;
        if sum(mnsSubGraphResults.breakingReaction(:,i),1) == 0 || isempty(listTemp)
            nlRange(i+1:end) = [];
            break;
        end
    end
end
% end
subGraphList = [num2cell([1:size(subGraphList,1)]') subGraphList];
if(initMNS.verboseScan == 2)
    disp('Scanning finished')
end
% idx = find(cell2mat(subGraphList(:,5)) > 1);
% tempList = subGraphList(idx,:);
%%
mnsSubGraphResults.sumBreakingReaction = sum(mnsSubGraphResults.breakingReaction,2);
mnsSubGraphResults.nlRange = nlRange;
mnsSubGraphResults.subGraphListHeader = {'SubGraphId','std0', 'std1', 'nL1' 'settingsSubGraphId', '#Edges', '#Nodes' '#Rp', '#Metabolites', 'Edges', 'EdgeIdx', ...
    'RpIdx', 'Model MetId', '#NoGenes', 'GeneIdx', '#Detected Metabolites','#unique Ions', 'Unique Ions Idx'};
mnsSubGraphResults.subGraphList = subGraphList;
 
% [mnsSubGraphResults.mergedSubGraphStruct, mnsSubGraphResults.mergedSubGraphList,...
%     mnsSubGraphResults.mergedSubGraphHeader] = merge_subGraphsByMet(mns, subGraphList, model.iaMat);
[mnsSubGraphResults.mergedSubGraphStruct, mnsSubGraphResults.mergedSubGraphList,...
    mnsSubGraphResults.mergedSubGraphHeader] = merge_subGraphs(mns, subGraphList, model.iaMat, initMNS.subGraphExtraction.mergeFreqCo);
%% calc p-value for all possible clusters...
if isfield(initMNS ,'determinePvalue') && initMNS.determinePvalue == 2
    if(initMNS.verboseScan == 2)
        disp('Start permutations: scan mode')
    end
    h = waitbar(0, ['0 of ' num2str(initMNS.permutations) ' permutations done']);
    initMnsTemp = initMNS;
    
    uSubGraphsTotalMet = unique(cell2mat(subGraphList(:,[1 3])), 'rows');
    uSubGraphsUniqueIons = unique(cell2mat(subGraphList(:,[1 6])), 'rows');
    freqSubGraphTotalMets = zeros(length(uSubGraphsTotalMet), initMnsTemp.permutations);
    freqSubGraphUniqueIons = zeros(length(uSubGraphsUniqueIons), initMnsTemp.permutations);
    isSubGraphTotalMets = zeros(length(uSubGraphsTotalMet), initMnsTemp.permutations);
    isSubGraphUniqueIons = zeros(length(uSubGraphsUniqueIons), initMnsTemp.permutations);
    
    modelTemp.iaMat = model.iaMat;
    modelTemp.metaboliteId = model.metaboliteId;
%     dataStructTemp = dataStruct;
    
    tp1 = toc(t1);
    if isfield(initMNS.parallel, 'parallelFolder') && ~isempty(initMNS.parallel.parallelFolder)
        mnsDataFolderP = initMNS.parallel.parallelFolder;
        mnsDataFolderCopy = initMNS.parallel.parallelFolderCopy;
        mnsExecFolderP = mnsDataFolderP;
    else
        mnsDataFolderP = mnsDataFolder;
        mnsDataFolderCopy = '';
        mnsExecFolderP = '';
    end
    if ~isempty(mnsDataFolderCopy)
        copyfile(fullfile(mnsFolder, 'opengm', mnsExec),...
            fullfile(mnsDataFolderCopy, mnsExec));
        copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'),...
            fullfile(mnsDataFolderCopy, 'nCliqueNoPerm.txt'));
        if(initMNS.input.edgeCliqueFileExists)
            copyfile(initMNS.input.edgeCliqueFile,...
                fullfile(mnsDataFolderCopy, 'cliques.txt'));
            initMnsTemp.input.edgeCliqueFile = fullfile(mnsDataFolderP, 'cliques.txt');
        end
    end
    copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'), fullfile(mnsDataFolderP, 'nCliqueNoPerm.txt'));
    if isfield(initMNS, 'parallel') && initMNS.parallel.useCluster
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
        elseif ~strcmp(poolobj.Cluster.Host, initMNS.parallel.clusterName) || initMNS.parallel.cores ~= poolobj.NumWorkers
            poolobj.delete();
            poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
        end
        stepSize = initMNS.parallel.cores*10;
        nLoops = ceil(initMNS.permutations/stepSize);
        sP = 0;
        pLoopSize = stepSize;
        for k = 1:nLoops
            const = (k-1)*stepSize;
            
            if sP+stepSize > initMnsTemp.permutations
                pLoopSize = initMnsTemp.permutations-sP;
            end
            parfor i = 1:pLoopSize
                dataStructTemp = dataStruct;
                dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
                mnsSubDataFolder = [mnsDataFolderP 'perm_' num2str(i)];
                [freqSubGraphTotalMets(:,const+i),freqSubGraphUniqueIons(:,const+i), ...
                    isSubGraphTotalMets(:,const+i), isSubGraphUniqueIons(:,const+i)] = ...
                    permStartMNS(modelTemp, initMnsTemp, dataStructTemp, ...
                    mnsDataFolderP, mnsSubDataFolder, nlRange,mnsExecFolderP,...
                    uSubGraphsTotalMet,uSubGraphsUniqueIons);
            end
            tp2 = toc(t1);
            sP = (k-1)*stepSize+pLoopSize;
%             if sP > initMnsTemp.permutation
%                sP = initMnsTemp.permutation;
%             end
            dt = (tp2-tp1)/sP;
            hours = floor(dt*(initMnsTemp.permutations-sP)/3600);
            mins = round((dt*(initMnsTemp.permutations-sP)-hours*3600)/60);
            waitbar(sP/initMNS.permutations, h, [{[num2str(sP) ' of ' num2str(initMNS.permutations) ' permutations done']}; {['Time remaining: ' num2str(hours) ' h ' num2str(mins) ' min ']}]);
        end
        
    else
        h = waitbar(0, ['0 of ' num2str(initMNS.permutations) ' permutations done']);
        for i = 1:initMnsTemp.permutations
            dataStructTemp = dataStruct;
            dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
            mnsSubDataFolder = [mnsDataFolderP 'perm_' num2str(i)];
            [freqSubGraphTotalMets(:,i),freqSubGraphUniqueIons(:,i), isSubGraphTotalMets(:,i), isSubGraphUniqueIons(:,i)] = ...
                permStartMNS(modelTemp, initMnsTemp, dataStructTemp, ...
                mnsDataFolderP, mnsSubDataFolder, nlRange,mnsExecFolderP,...
                uSubGraphsTotalMet,uSubGraphsUniqueIons);
            tp2 = toc(t1);
            dt = (tp2-tp1)/i;
            hours = floor(dt*(initMnsTemp.permutations-i)/3600);
            mins = round((dt*(initMnsTemp.permutations-i)-hours*3600)/60);
            waitbar(i/initMNS.permutations, h, [{[num2str(i) ' of ' num2str(initMNS.permutations) ' permutations done']}; {['Time remaining: ' num2str(hours) ' h ' num2str(mins) ' min ']}]);
        end
        
    end
    delete(h);
    
    
    mnsSubGraphResults.p.pMode = 'scan mode';
    mnsSubGraphResults.p.combSubGraphsTotalMetHeader = {'nL1' '#Metabolites'};
    mnsSubGraphResults.p.combSubGraphsTotalMet = uSubGraphsTotalMet;
    mnsSubGraphResults.p.pSubGraphTotalMets = sum(isSubGraphTotalMets,2)/size(isSubGraphTotalMets,2);
    mnsSubGraphResults.p.freqSubGraphTotalMets = sum(freqSubGraphTotalMets,2);
    
    mnsSubGraphResults.p.combSubGraphsUniqueIonsHeader = {'nL1' '#Unique Ions'};
    mnsSubGraphResults.p.combSubGraphsUniqueIons = uSubGraphsUniqueIons;
    mnsSubGraphResults.p.pSubGraphUniqueIons = sum(isSubGraphUniqueIons,2)/size(isSubGraphUniqueIons,2);
    mnsSubGraphResults.p.freqSubGraphUniqueIons = sum(freqSubGraphUniqueIons,2);
    if(initMNS.verboseScan == 2)
        disp('Permutations done')
    end
end



try
    rmdir(mnsDataFolder, 's');
catch err
    disp(err)
end
t2 = toc(t1);
if(initMNS.verboseScan ~= 0)
    disp(['elapsed time: ' num2str(floor(t2/60)) ' min ' num2str(round(mod(t2,60))) ' sec']);
end
end

%% function to get all subgraphs given the breaking reactions....
function [subGraphs, subGraphList] = getSubgraphs(breakingReaction, ...
    reactionPairs, iaMat, iaId2DataId, inferedLabels, isRp)
rpReduced = reactionPairs(breakingReaction == 1,:);
for i = 1:size(rpReduced,1)
    iaMat(rpReduced(i,1), rpReduced(i,2)) = 0;
    iaMat(rpReduced(i,2), rpReduced(i,1)) = 0;
end

sPaths = graphallshortestpaths(sparse(iaMat), 'Directed', false);
subGraphs.cmpdSubGraphLabels = ones(size(iaMat,1),1)*-1;
cMetLabel = 0;
subGraphList = [];
while(min(subGraphs.cmpdSubGraphLabels) == -1)
    idx = find(subGraphs.cmpdSubGraphLabels == -1);
    idx = idx(1);
    if inferedLabels(idx) == 0 
        subGraphs.cmpdSubGraphLabels(isinf(sPaths(idx,:)) == 0) = 0;
    else
        tempMetIdx = find(isinf(sPaths(idx,:)) == 0);
        tempNoOfMet = 0;
        tempNoOfDetectedMet = 0;
        tempNoOfRp = 0;
        tempMet = [];
        tempDetectedMet = [];
        tempGene = [];
        for i = 1:length(tempMetIdx)
            if isRp(tempMetIdx(i))
                temp = iaId2DataId(tempMetIdx(i));
                tempNoOfRp = tempNoOfRp+1;
                if ~isempty(temp);
                    
                    tempGene = [tempGene; temp];
                end
            else
                temp = iaId2DataId{tempMetIdx(i)};
                tempNoOfMet = tempNoOfMet+1;
                tempMet = [tempMet; tempMetIdx(i)];
                if ~isempty(temp)
                    tempNoOfDetectedMet = tempNoOfDetectedMet+1;
                    tempDetectedMet = [tempDetectedMet; temp];
                end
            end
        end
        if ~isempty(tempDetectedMet) && tempNoOfDetectedMet >= 2
            cMetLabel = cMetLabel+1;
            
            %get model indices for cmpd (rp+met) met and rp
            subGraphs.subGraphDetails(cMetLabel).cmpdIdx = find(isinf(sPaths(idx,:)) == 0);
            subGraphs.subGraphDetails(cMetLabel).metIdx = tempMet;
            subGraphs.subGraphDetails(cMetLabel).rpIdx = setdiff(subGraphs.subGraphDetails(cMetLabel).cmpdIdx,tempMet);
            
            % get edges
            edges = [];
            edgeIdx = [];
            for i = 1:length(subGraphs.subGraphDetails(cMetLabel).cmpdIdx)              
                idxTemp = find(iaMat(subGraphs.subGraphDetails(cMetLabel).cmpdIdx(i),:) ~= 0);
                idxTemp(idxTemp <= subGraphs.subGraphDetails(cMetLabel).cmpdIdx(i)) = [];
                edges = [edges; ones(length(idxTemp),1)*subGraphs.subGraphDetails(cMetLabel).cmpdIdx(i) idxTemp'];
            end
            edgeIdx = [(edges(:,1)-1).*size(iaMat,1)+edges(:,2)];
            
            %get sizes
            subGraphs.subGraphSize(cMetLabel) = length(subGraphs.subGraphDetails(cMetLabel).cmpdIdx);
            subGraphs.subGraphSizeMet(cMetLabel) = length(subGraphs.subGraphDetails(cMetLabel).metIdx);
            subGraphs.subGraphSizeRp(cMetLabel) = length(subGraphs.subGraphDetails(cMetLabel).rpIdx);
            subGraphs.subGraphSizeEdges(cMetLabel) = length(edgeIdx);
            subGraphs.subGraphSizeGenes(cMetLabel) = length(tempGene);
            
            %edge information
            subGraphs.subGraphDetails(cMetLabel).edges = edges;
            subGraphs.subGraphDetails(cMetLabel).edgeIdx = edgeIdx;
            
            %cmpd information
            subGraphs.cmpdSubGraphLabels(subGraphs.subGraphDetails(cMetLabel).cmpdIdx) = cMetLabel;
            
            %get unique ions
            subGraphs.subGraphDetails(cMetLabel).dataGeneIdx = tempGene;
            subGraphs.subGraphDetails(cMetLabel).dataMetIdx = unique(tempDetectedMet);
            subGraphs.subGraphUniqueIons(cMetLabel) = length(subGraphs.subGraphDetails(cMetLabel).dataMetIdx);
            subGraphs.subGraphDetectedMetabolites(cMetLabel) = tempNoOfDetectedMet;
            
            subGraphList = [subGraphList; cMetLabel subGraphs.subGraphSizeEdges(cMetLabel)...
                subGraphs.subGraphSize(cMetLabel) subGraphs.subGraphSizeRp(cMetLabel) ...
                subGraphs.subGraphSizeMet(cMetLabel) ...
                {subGraphs.subGraphDetails(cMetLabel).edges} {subGraphs.subGraphDetails(cMetLabel).edgeIdx}...
                {subGraphs.subGraphDetails(cMetLabel).rpIdx} {subGraphs.subGraphDetails(cMetLabel).metIdx} ...
                subGraphs.subGraphSizeGenes(cMetLabel) {subGraphs.subGraphDetails(cMetLabel).dataGeneIdx} ...
                subGraphs.subGraphDetectedMetabolites(cMetLabel) subGraphs.subGraphUniqueIons(cMetLabel)...
                {subGraphs.subGraphDetails(cMetLabel).dataMetIdx } ...
                ];
        else
            subGraphs.cmpdSubGraphLabels(isinf(sPaths(idx,:)) == 0) = 0;
        end
    end
end
% subGraphs.cmpdSubGraphLabels((subGraphs.cmpdSubGraphLabels) == -1) = 0;

end

%% merging of subGraphs to find overlapping subgraphs

function [mergedSubGraphStruct, mergedSubGraphList, mergedSubGraphHeader] ...
    = merge_subGraphs(mns, subGraphList, iaMat,freqCo)
%% predefine matrices and get data
if isempty(subGraphList)
    mergedSubGraphStruct = [];
    mergedSubGraphList = [];
    mergedSubGraphHeader = [];
    return;
end

intersectMat = zeros(size(subGraphList,1)); %overlap of edges in subgraphs
intersectMatEdgeList = cell(size(subGraphList,1));
mergedSubGraphStruct.cmpdSubGraphLabels = zeros(length(mns.observations.modelVarId),1);
noEdges = cell2mat(subGraphList(:,6));
edgeIdx = subGraphList(:,11);
nCmpd = size(iaMat,1);
%% check for overlap between edges

for i = 1:size(subGraphList,1)
    intersectMat(i,i) = noEdges(i);
    for j = i+1:size(subGraphList,1)
        intersectMatEdgeList(i,j) = {intersect(edgeIdx{i},edgeIdx{j})};
        if ~isempty(intersectMatEdgeList(i,j))
            intersectMat(i,j) = length(intersectMatEdgeList{i,j});
            intersectMat(j,i) = length(intersectMatEdgeList{i,j});
        end
    end
end

%% make iaMat and find connected subgraphs
iaMatEdges = zeros(size(intersectMat));
iaMatEdges(intersectMat ~= 0) = 1;


%% merge graphs for edges
sPaths = graphallshortestpaths(sparse(iaMatEdges), 'Directed', false);
subGraphLabels = zeros(length(edgeIdx),1);
mergedSubGraphList = [];
c = 0;
%%
while(nanmin(subGraphLabels) == 0)
    clear tempStruct
    idx = find(subGraphLabels == 0);
    idx = idx(1);
    tempStruct.tempGraphIdx = find(isinf(sPaths(idx,:)) == 0);
    subGraphLabels(tempStruct.tempGraphIdx) = nan(1);
    [tempStruct.edgeIdxUnion, tempStruct.edgeIdxFrequency] = ...
        subGraphUnion(edgeIdx(tempStruct.tempGraphIdx));
    if min(tempStruct.edgeIdxFrequency) < freqCo
        tempGraphIdx = tempStruct.tempGraphIdx;
        tempStruct = splitMergedSubGraphsByCo(subGraphList(tempGraphIdx,:), ...
            tempStruct.edgeIdxUnion, tempStruct.edgeIdxFrequency, freqCo,size(iaMat,1));
        if isempty(tempStruct)
            continue;
        end
        for i = 1:length(tempStruct)
            tempStruct(i).tempGraphIdx(tempStruct(i).tempGraphIdx > length(tempGraphIdx)) = [];
%             disp(i)
            tempStruct(i).tempGraphIdx = tempGraphIdx(tempStruct(i).tempGraphIdx);
        end
    end
    %%
    for i = 1:length(tempStruct)
        if isempty(tempStruct(i).edgeIdxFrequency)
            continue;
        end
        tempGraphIdx = tempStruct(i).tempGraphIdx;
        edgeIdxFrequency = tempStruct(i).edgeIdxFrequency;
        edgeIdxUnion = tempStruct(i).edgeIdxUnion;
        c = c+1;
        subGraphLabels(tempGraphIdx) = c;
        
        %store edge information
        mergedSubGraphStruct.subGraphEdgeMeanFreq(c)= mean(edgeIdxFrequency);
        mergedSubGraphStruct.subGraphDetails(c).edgeIdx = edgeIdxUnion;
        mergedSubGraphStruct.subGraphDetails(c).edgeIdxFrequency = edgeIdxFrequency;
        mergedSubGraphStruct.subGraphDetails(c).edgeIdxCore = edgeIdxUnion(edgeIdxFrequency == 1);
        mergedSubGraphStruct.subGraphDetails(c).edges = [(edgeIdxUnion-mod(edgeIdxUnion,nCmpd))/nCmpd+1, mod(edgeIdxUnion,nCmpd)];
        idxEdgeTemp = find(mergedSubGraphStruct.subGraphDetails(c).edges(:,2) == 0);
        if ~isempty(idxEdgeTemp)
            mergedSubGraphStruct.subGraphDetails(c).edges(idxEdgeTemp,1) = mergedSubGraphStruct.subGraphDetails(c).edges(idxEdgeTemp,1)-1;
            mergedSubGraphStruct.subGraphDetails(c).edges(idxEdgeTemp,2) = nCmpd;
        end
        mergedSubGraphStruct.subGraphSizeEdges(c) = length(edgeIdxUnion);
        mergedSubGraphStruct.coreSubGraphSize(c) = length(mergedSubGraphStruct.subGraphDetails(c).edgeIdxCore);
        mergedSubGraphStruct.averageEdgeFrequency(c) = mean(edgeIdxFrequency);
        
        % get cmpdIdx from edges
        mergedSubGraphStruct.subGraphDetails(c).cmpdIdx = ...
            unique(reshape(mergedSubGraphStruct.subGraphDetails(c).edges,length(edgeIdxUnion)*2,1));
        geneDataIdx = [];
        metDataIdx = [];
        metIdx = [];
        rpIdx = [];
        geneIdx = [];
        for j = 1:length(mergedSubGraphStruct.subGraphDetails(c).cmpdIdx)
            if mergedSubGraphStruct.subGraphDetails(c).cmpdIdx(j) > length(mns.observations.isRp)
                mergedSubGraphStruct.subGraphDetails(c).cmpdIdx
            end
            if mns.observations.isRp(mergedSubGraphStruct.subGraphDetails(c).cmpdIdx(j))
                rpIdx = [rpIdx; mergedSubGraphStruct.subGraphDetails(c).cmpdIdx(j)];
                tempAnnotation = mns.observations.iaId2dataId(mergedSubGraphStruct.subGraphDetails(c).cmpdIdx(j));
                if ~isempty(tempAnnotation)
                    geneIdx = [geneIdx; tempAnnotation];
                end
            else
                metIdx = [metIdx; mergedSubGraphStruct.subGraphDetails(c).cmpdIdx(j)];
                tempAnnotation = mns.observations.iaId2dataId{mergedSubGraphStruct.subGraphDetails(c).cmpdIdx(j)};
                if tempAnnotation ~= 0
                    metDataIdx = [metDataIdx; tempAnnotation];
                end
            end 
        end
        
        %get met, rp and geneIdx
        mergedSubGraphStruct.subGraphDetails(c).metIdx = metIdx;
        mergedSubGraphStruct.subGraphDetails(c).rpIdx = setdiff(mergedSubGraphStruct.subGraphDetails(c).cmpdIdx,metIdx);

        %define sizes
        mergedSubGraphStruct.subGraphSize(c) = length(mergedSubGraphStruct.subGraphDetails(c).cmpdIdx);
        mergedSubGraphStruct.subGraphSizeMet(c) = length(mergedSubGraphStruct.subGraphDetails(c).metIdx);
        mergedSubGraphStruct.subGraphSizeRp(c) = length(mergedSubGraphStruct.subGraphDetails(c).rpIdx);
%         subGraphs.subGraphSizeEdges(c) = length(edgeIdx);
        mergedSubGraphStruct.subGraphSizeGenes(c) = length(geneIdx);
        
        %get detected metabolites unique ions
        mergedSubGraphStruct.subGraphDetails(c).geneIdx = geneIdx;
        mergedSubGraphStruct.subGraphDetails(c).dataMetIdx = unique(metDataIdx);
        mergedSubGraphStruct.subGraphUniqueIons(c) = length(mergedSubGraphStruct.subGraphDetails(c).dataMetIdx);
        mergedSubGraphStruct.subGraphDetectedMetabolites(c) = length(metDataIdx);

        mergedSubGraphList = [mergedSubGraphList; ...
            c mergedSubGraphStruct.subGraphSizeEdges(c) mergedSubGraphStruct.coreSubGraphSize(c) ...
            mergedSubGraphStruct.subGraphSize(c) mergedSubGraphStruct.subGraphSizeRp(c) ...
            mergedSubGraphStruct.subGraphSizeMet(c) ...
            mean(edgeIdxFrequency) {mergedSubGraphStruct.subGraphDetails(c).edges} ...
            {edgeIdxUnion} {edgeIdxFrequency} ...
            {mergedSubGraphStruct.subGraphDetails(c).rpIdx} {mergedSubGraphStruct.subGraphDetails(c).metIdx}...
            mergedSubGraphStruct.subGraphSizeGenes(c) {mergedSubGraphStruct.subGraphDetails(c).geneIdx}...
            mergedSubGraphStruct.subGraphDetectedMetabolites(c) mergedSubGraphStruct.subGraphUniqueIons(c)...$
            {mergedSubGraphStruct.subGraphDetails(c).dataMetIdx}
            ];
    end
end
 mergedSubGraphHeader = {'SubGraphIdx' '#Edges' '#Core Edges' '#Nodes' '#Rp', '#Metabolites', 'AvgFreq' ...
     'Edges' 'EdgeIdx' 'EdgeFreq' 'RpIdx', 'MetIdx', '#NoGenes', ...
     'GeneIdx', '#Detected Metabolites','#unique Ions', 'Unique Ions Idx'...
     };
end

%% merging of subGraphs by metabolites to find overlapping subgraphs

function [mergedSubGraphStruct, mergedSubGraphList, mergedSubGraphHeader] ...
    = merge_subGraphsByMet(mns, subGraphList, iaMat)
%% predefine matrices and get data
intersectMat = zeros(size(subGraphList,1)); %overlap of metabolites in subgraphs
intersectMatMetList = cell(size(subGraphList,1));

mergedSubGraphStruct.cmpdSubGraphLabels = zeros(length(mns.observations.metId),1);
noMet = cell2mat(subGraphList(:,6));
metIdx = subGraphList(:,7);
freqCo = 0.2;
%% check for overlap between metabolites
for i = 1:length(subGraphList)
    intersectMat(i,i) = noMet(i);
%     intersectMatUniqueIons(i,i) = noUniqueIons(i);
    for j = i+1:length(subGraphList)
        intersectMatMetList(i,j) = {intersect(metIdx{i},metIdx{j})};
%         intersectMatUniqueIonsList{i,j} = intersect(ionIdx{i},ionIdx{j});
        if ~isempty(intersectMatMetList(i,j))
            intersectMat(i,j) = length(intersectMatMetList{i,j})/max(length(metIdx{i}),length(metIdx{j}));
            intersectMat(j,i) = length(intersectMatMetList{i,j})/max(length(metIdx{i}),length(metIdx{j}));
        end
%         if ~isempty(intersectMatUniqueIonsList(i,j))
%             intersectMatUniqueIons(i,j) = length(intersectMatUniqueIonsList{i,j});
%             intersectMatUniqueIons(j,i) = length(intersectMatUniqueIonsList{i,j});
%         end
    end
end
%% make iaMat and find connected subgraphs
iaMatMet = zeros(size(intersectMat));
iaMatMet(intersectMat > freqCo) = 1;
sPaths = graphallshortestpaths(sparse(iaMatMet), 'Directed', false);

%% merge graphs
subGraphLabels = zeros(length(metIdx),1);
c = 0;
mergedSubGraphList = [];
while(min(subGraphLabels) == 0)
    c = c+1;
    idx = find(subGraphLabels == 0);
    idx = idx(1);
    tempGraphIdx = find(isinf(sPaths(idx,:)) == 0);
    subGraphLabels(tempGraphIdx) = c;
    [metIdxUnion, metIdxFrequency] = subGraphUnionByMet(metIdx(tempGraphIdx));
    mergedSubGraphStruct.cmpdSubGraphLabels(metIdxUnion) = c;
    mergedSubGraphStruct.subGraphDetails(c).metIdx = metIdxUnion;
    mergedSubGraphStruct.subGraphDetails(c).metIdxFrequency = metIdxFrequency;
    mergedSubGraphStruct.subGraphDetails(c).metIdxCore = metIdxUnion(metIdxFrequency == 1);
    mergedSubGraphStruct.subGraphSize(c) = length(metIdxUnion);
    mergedSubGraphStruct.coreSubGraphSize(c) = length(mergedSubGraphStruct.subGraphDetails(c).metIdxCore);
    mergedSubGraphStruct.averageMetFrequency(c) = mean(metIdxFrequency);
    ionIdxUnion = mns.observations.iaId2dataId(metIdxUnion);
    mergedSubGraphStruct.subGraphDetails(c).modelMetIdx2dataIdx = ionIdxUnion;
        
    ionIdxUnion(ionIdxUnion == 0) = 0;
    mergedSubGraphStruct.subGraphDetectedMetabolites(c) = length(ionIdxUnion);
    mergedSubGraphStruct.subGraphDetails(c).ionIdxUnion = unique(ionIdxUnion);
    mergedSubGraphStruct.subGraphUniqueIons(c) = length(mergedSubGraphStruct.subGraphDetails(c).ionIdxUnion);
    
    mergedSubGraphList = [mergedSubGraphList; ...
        num2cell([c mergedSubGraphStruct.subGraphSize(c) ...
        mergedSubGraphStruct.coreSubGraphSize(c) mean(metIdxFrequency)]) ...
        {metIdxUnion} {metIdxFrequency} ...
        num2cell([mergedSubGraphStruct.subGraphDetectedMetabolites(c)...
        mergedSubGraphStruct.subGraphUniqueIons(c)])...
        {mergedSubGraphStruct.subGraphDetails(c).ionIdxUnion}];
end
mergedSubGraphHeader = {'SubGraphIdx' '#Metabolites' '#Core Metabolites' 'AvgFreq' 'MetIdx'...
    'MetFreq' '#Detected Metabolites', '#unique Ions' 'ionIdx'};
end

%% make unione of subgraph metabolites
function [metIdxUnion, metIdxFrequency] = subGraphUnionByMet(metStruct)
    metIdxUnion = metStruct{1};
    metIdxArr = metStruct{1};
    for i = 2:length(metStruct)
        metIdxUnion = union(metIdxUnion, metStruct{i});
        metIdxArr = [metIdxArr metStruct{i}];
    end
    if length(metIdxUnion) == 1
        [metIdxFrequency,b] = hist(metIdxArr, 1);
    else
        [metIdxFrequency,b] = hist(metIdxArr, metIdxUnion);
    end
    
    metIdxFrequency = metIdxFrequency/length(metStruct);
end

%% make unione of subgraph metabolites
function [edgeIdxUnion, edgeIdxFrequency] = subGraphUnion(edgeStruct)
    edgeIdxUnion = edgeStruct{1};
    edgeIdxArr = edgeStruct{1}';
    for i = 2:length(edgeStruct)
        edgeIdxUnion = union(edgeIdxUnion, edgeStruct{i});
        edgeIdxArr = [edgeIdxArr edgeStruct{i}'];
    end
    if length(edgeIdxUnion) == 1
        [edgeIdxFrequency,b] = hist(edgeIdxArr, 1);
    else
        [edgeIdxFrequency,b] = hist(edgeIdxArr, edgeIdxUnion);
    end
    
    edgeIdxFrequency = edgeIdxFrequency/length(edgeStruct);
end

function outStruct = splitMergedSubGraphsByCo(subGraphList, edgeIdxUnion, edgeIdxFrequency, freqCo, noOfMet)
%% predefine matrices and get data
noEdges = cell2mat(subGraphList(:,6));
edgeIdx = subGraphList(:,11);
edges = subGraphList(:,10);

%% check for overlap between edges
edgesToRemove = edgeIdxUnion(edgeIdxFrequency < freqCo);
for i = 1:size(subGraphList,1)
    [~,idxA] = intersect(edgeIdx{i}, edgesToRemove);
    edgeIdx{i}(idxA) = [];
    edges{i}(idxA,:) = [];
    %% check if subgraph still exists
    iaMatTemp = zeros(noOfMet);
    for j = 1:size(edges{i},1)
%         iaMatTemp(edges{i}(j,1),edges{i}(j,2)) = 1;
        iaMatTemp(edges{i}(j,2),edges{i}(j,1)) = 1;
    end
    uMet = unique([edges{i}(:,1);edges{i}(:,2)]);
    sPathsTemp = graphallshortestpaths(sparse(iaMatTemp), 'Directed', false);
    sPathsTemp = sPathsTemp(uMet,uMet);
%     iaMatTemp = iaMatTemp(uMet,uMet);
    %%
    firstCluster = true;
    firstEdge = true;
    clLabel = zeros(length(uMet),1);
    %%
    while(min(clLabel) == 0)
        %%
        idxTempMet = find(clLabel == 0);
        j = idxTempMet(1);
        idxTemp = find(isinf(sPathsTemp(j,:)) == 0);
        if firstCluster
            c = i;
            firstCluster = false;
        else
            c = size(edges,1)+1;
            firstEdge = true;
        end
        clLabel(idxTemp) = c;
        idxTemp = uMet(idxTemp);
        if ~isempty(idxTemp)
            for k = 1:length(idxTemp)
                for l = k+1:length(idxTemp)
                    if iaMatTemp(idxTemp(k),idxTemp(l)) || iaMatTemp(idxTemp(l),idxTemp(k))
                        minEdge = min([idxTemp(k) idxTemp(l)]);
                        maxEdge = max([idxTemp(k) idxTemp(l)]);
                        if firstEdge
                            edges{c} = [minEdge maxEdge];
                            edgeIdx{c} = (minEdge-1)*noOfMet + maxEdge;
                            noEdges(c) = 1;
                            firstEdge = false;
                            
                        else
                            if ~isempty(find(edges{c}(:,1) == minEdge & edges{c}(:,2) == maxEdge))
                                continue;
                            end
                            edges{c} = [edges{c}; minEdge maxEdge];
                            edgeIdx{c} = [edgeIdx{c};(minEdge-1)*noOfMet + maxEdge];
                            noEdges(c) = noEdges(c)+1;
                        end
                    end
                end
            end
            
        end
    end
end

intersectMat = zeros(length(edges)); %overlap of edges in subgraphs
intersectMatEdgeList = cell(length(edges));

%%
for i = 1:length(edgeIdx)
    intersectMat(i,i) = noEdges(i);
    for j = i+1:length(edgeIdx)
        intersectMatEdgeList(i,j) = {intersect(edgeIdx{i},edgeIdx{j})};
        if ~isempty(intersectMatEdgeList(i,j))
            intersectMat(i,j) = length(intersectMatEdgeList{i,j});
            intersectMat(j,i) = length(intersectMatEdgeList{i,j});
        end
    end
end

%% make iaMat and find connected subgraphs
iaMatEdges = zeros(size(intersectMat));
iaMatEdges(intersectMat ~= 0) = 1;
sPaths = graphallshortestpaths(sparse(iaMatEdges), 'Directed', false);

%% merge graphs for edges
% mergedSubGraphList = [];
subGraphLabels = zeros(length(edgeIdx),1);
c = 0;
outStruct = [];
while(nanmin(subGraphLabels) == 0)
    
    idx = find(subGraphLabels == 0);
    idx = idx(1);
    if length(find(isinf(sPaths(idx,:)) == 0)) > 1
        c = c+1;
        outStruct(c).tempGraphIdx = find(isinf(sPaths(idx,:)) == 0);
        subGraphLabels(outStruct(c).tempGraphIdx) = c;
        [outStruct(c).edgeIdxUnion, outStruct(c).edgeIdxFrequency] = ...
            subGraphUnion(edgeIdx(outStruct(c).tempGraphIdx));
    else
        subGraphLabels(idx) = nan(1);
    end
end

end
%% function for p-value determiniation (needs to be adapted to std scanning)
function [freqSubGraphTotalMets, freqSubGraphUniqueIons, isSubGraphTotalMets, isSubGraphUniqueIons]...
    = permStartMNS(model, initMNS, dataStruct, mnsFolder, mnsSubDataFolderStem, nlRange,mnsExecFolderP,uSubGraphsTotalMet,uSubGraphsUniqueIons)
    
    subGraphList = [];
    for i = 1:length(nlRange)
        mnsSubDataFolder = [mnsSubDataFolderStem '_nl_' num2str(nlRange(i)) '\'];
        mkdir(mnsSubDataFolder);
        copyfile([mnsFolder 'nCliqueNoPerm.txt'], [mnsSubDataFolder 'nCliqueNoPerm.txt']);
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, [mnsSubDataFolder 'cliques.txt']);
        end
        [mnsTemp, ~, ~] = mns_run2state(model, dataStruct, initMNS, mnsSubDataFolder,mnsExecFolderP);
        [~, tempList] = getSubgraphs(mnsTemp.results.breakingReactions, mnsTemp.results.reactionPairs, ...
            model.iaMat, mnsTemp.observations.iaId2dataId, mnsTemp.results.inferedLabels);
        subGraphList = [subGraphList; num2cell(ones(size(tempList,1),1).*nlRange(i)) tempList];
    end

    freqSubGraphTotalMets = zeros(length(uSubGraphsTotalMet), 1);
    freqSubGraphUniqueIons = zeros(length(uSubGraphsUniqueIons), 1);
    isSubGraphTotalMets = zeros(length(uSubGraphsTotalMet), 1);
    isSubGraphUniqueIons = zeros(length(uSubGraphsUniqueIons), 1);
    if ~isempty(subGraphList)
        for i = 1:size(uSubGraphsTotalMet,1)
            idx = find(cell2mat(subGraphList(:,1)) == uSubGraphsTotalMet(i,1) ...
                & cell2mat(subGraphList(:,3)) >= uSubGraphsTotalMet(i,2));
            if ~isempty(idx)
                freqSubGraphTotalMets(i) = length(idx);
                isSubGraphTotalMets(i) = length(idx);
            end
        end
        for i = 1:size(uSubGraphsUniqueIons,1)
            idx = find(cell2mat(subGraphList(:,1)) == uSubGraphsUniqueIons(i,1) ...
                & cell2mat(subGraphList(:,6)) >= uSubGraphsUniqueIons(i,2));
            if ~isempty(idx)
                freqSubGraphUniqueIons(i) = length(idx);
                isSubGraphUniqueIons(i) = length(idx);
            end
        end
        isSubGraphUniqueIons(isSubGraphUniqueIons ~= 0) = 1;
        isSubGraphTotalMets(isSubGraphTotalMets ~= 0) = 1;
    end
    try
        rmdir(mnsSubDataFolder,'s');
    catch err
        disp(err)
    end
end