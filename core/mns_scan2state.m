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
function mnsScanResults = mns_scan2state(model,dataStruct, nameTag, nlRange, initMNS, plotResults)
t1 = tic;
if nargin < 3 || isempty(nameTag)
    nameTag = 'temp';
end
if nargin < 4 || isempty(nlRange)
    nlRange = 0;
end
if nargin < 5 || isempty(initMNS)
    initMNS = mns_generateInitMNS;
end
if nargin < 6 || isempty(plotResults)
    plotResults = false;
end
if isfield(dataStruct, 'dataReplicates')
    initMNS.biolReplicates = dataStruct.dataReplicates;
end
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
    otherwise
end
if ~isempty(initMNS.input.dataFolder) && isfield(initMNS.input, 'dataFolder')
    mnsDataFolderStem = initMNS.input.dataFolder;
    clFolder = mnsDataFolderStem;
%     copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'),fullfile(clFolder, 'nCliqueNoPerm.txt'));
    if ~isempty(initMNS.input.mnsExecFolder) && isfield(initMNS.input, 'mnsExecFolder')
        
%         copyfile(fullfile(mnsFolder, 'opengm', mnsExec),...
%             fullfile(initMNS.input.folderCopy, mnsExec));
%         copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'),...
%             fullfile(initMNS.input.folderCopy, 'nCliqueNoPerm.txt'));
        mnsExecFolder = initMNS.input.mnsExecFolder;
    else
        mnsExecFolder = '';
    end
else
    clFolder = mnsFolder;
    mnsDataFolderStem = fullfile(mnsBaseFolder, 'data');
%     switch computer
%         case 'PCWIN32'
%             mnsDataFolderStem = [mnsFolder '\data'];
%         case 'PCWIN64'
%             mnsDataFolderStem = [mnsFolder '\data'];
%         otherwise
%             mnsDataFolderStem = [mnsFolder '/data'];
%     end
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
mnsDataFolder = fullfile(mnsDataFolderStem,nameTag);
% disp(mnsDataFolderStem)
% disp(nameTag)
switch computer
    case 'PCWIN32'
        mnsDataFolder = [mnsDataFolder '\'];
%         regexprep(mnsDataFolder,'\\','\\\\');
    case 'PCWIN64'
        mnsDataFolder = [mnsDataFolder '\'];
%         regexprep(mnsDataFolder,'\\','\\\\');
    otherwise
        mnsDataFolder = [mnsDataFolder '/'];
end

mkdir(mnsDataFolder)

% if nlRange = 0, do coarse grained search
if nlRange == 0
    if(initMNS.verboseScan == 2)
        disp('Start coarse grained search...');
    end
    nlTemp = 1;
    initMNStemp = initMNS;
    i = 0;
    while(nlTemp <= 512)
        i = i+1;
        mnsSubDataFolderTemp = fullfile(mnsDataFolder, ['coarse_nl_' num2str(nlTemp)]);
        mkdir(mnsSubDataFolderTemp);
        copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(mnsSubDataFolderTemp, 'nCliqueNoPerm.txt'));
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, fullfile(mnsSubDataFolderTemp, 'cliques.txt'));
        end
        if(initMNS.verboseScan == 2)
            disp(nlTemp);
        end
        initMNStemp.nL1max = nlTemp;
        initMNStemp.nL1min= nlTemp;
        
        [mns] = mns_run2state(model, dataStruct, initMNStemp, mnsSubDataFolderTemp, mnsExecFolder);
        if isempty(mns.results)
            mnsScanResults = [];
            return;
        end
        if i == 1
            initMNStemp.input.initLabels = mns.results.initLabels;
            initMNStemp.labeling = 'fix';
            initMNStemp.input.reactionPairs = mns.results.reactionPairs;
            initMNStemp.input.edgeCliqueFile = fullfile(mnsDataFolder,  ['coarse_nl_' num2str(nlTemp)], 'cliques.txt');
            initMNStemp.input.edgeCliqueFileExists = true;
            initMNStemp.input.nCliques = mns.nCliques;
            initMNS.input.edgeCliqueFile = fullfile(mnsDataFolder,  ['coarse_nl_' num2str(nlTemp)], 'cliques.txt');
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
fullNlRange = true;
for i = 1:length(nlRange)
    mnsSubDataFolderTemp = fullfile(mnsDataFolder, ['nl_' num2str(nlRange(i))]);
    mkdir(mnsSubDataFolderTemp);
    copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(mnsSubDataFolderTemp, 'nCliqueNoPerm.txt'));
    if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
        copyfile(initMNS.input.edgeCliqueFile, fullfile(mnsSubDataFolderTemp,'cliques.txt'));
    end
    if(initMNS.verboseScan == 2)
        disp('');
        disp(['Starting step ' num2str(i) ', nl = ' num2str(nlRange(i))]);
        disp('');
    end
    initMNS.nL1max = nlRange(i);
    initMNS.nL1min= nlRange(i);
    [mns] = mns_run2state(model, dataStruct, initMNS, mnsSubDataFolderTemp, mnsExecFolder);
    if isempty(mns.results)
        mnsScanResults = [];
        return;
    end
    mnsScanResults.mns(i) = {mns};
    mnsScanResults.breakingReaction(:,i) = mns.results.breakingReactions;
    mnsScanResults.clustLabel(:,i) = mns.results.inferedLabels;
    
    if i == 1
        initMNS.input.edgeCliqueFile = fullfile(mnsDataFolder, ['nl_' num2str(nlRange(i))], 'cliques.txt');
        initMNS.input.edgeCliqueFileExists = true;
        initMNS.input.nCliques = mns.nCliques;
        initMNS.input.initLabels = mnsScanResults.mns{1}.results.initLabels;
        initMNS.labeling = 'fix';
        initMNS.input.reactionPairs = mnsScanResults.mns{1}.results.reactionPairs;
        mnsScanResults.maxBreakingReaction = zeros(size(mnsScanResults.breakingReaction,1),1);
        mnsScanResults.sumBreakingReaction = zeros(size(mnsScanResults.breakingReaction,1),1);
    end
    mnsScanResults.maxBreakingReaction(mnsScanResults.breakingReaction(:,i) == 1) = i;
    if sum(mnsScanResults.breakingReaction(:,i),1) == 0
        nlRange(i+1:end) = [];
        fullNlRange = false;
        break;
    end
end
% end
if(initMNS.verboseScan == 2)
    disp('Scanning finished')
end
mnsScanResults.absDiffSubProdReaction = getAbsDiffSubProdReaction(mnsScanResults.mns{1}.results.reactionPairs,mnsScanResults.mns{1}.observations.data);
mnsScanResults.sumBreakingReaction = sum(mnsScanResults.breakingReaction,2);
mnsScanResults.maxSumBreakingReaction = mnsScanResults.sumBreakingReaction.*mnsScanResults.maxBreakingReaction;
mnsScanResults.nlRange = nlRange;
if fullNlRange
    nlRange(end+1) = nan;
end
%% calc p-value only for maxRecations
if isfield(initMNS ,'determinePvalue') && initMNS.determinePvalue == 1
    if(initMNS.verboseScan == 2)
        disp('Start permutations: only max')
    end
    nl = nlRange(end-1);
    initMnsTemp = initMNS;
    initMnsTemp.nL1max = nlRange(end-1);
    initMnsTemp.nL1min = nlRange(end-1);
    noOfBreakingPermutations = zeros(initMNS.permutations,1);
    modelTemp.iaMat = model.iaMat;
    modelTemp.metaboliteId = model.metaboliteId;
    
    if isfield(initMNS.parallel, 'useCluster') && initMNS.parallel.useCluster
        poolobj = gcp('nocreate');
        if isempty(poolobj) || ~strcmp(initMNS.parallel.clusterName, poolobj.Cluster.profile) || initMNS.parallel.cores ~= poolobj.Cluster.NumWorkers
            poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
        end
        parfor i = 1:initMnsTemp.permutations
            dataStructTemp = dataStruct;
            dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
            mnsSubDataFolder = [mnsDataFolder 'perm_' num2str(i) '_nl_' num2str(nl) '\'];
            [breakingReactionsPerm(:,i), noOfBreakingPermutations(i)] = ...
                permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsFolder, mnsSubDataFolder);
        end
    else
        for i = 1:initMnsTemp.permutations
            dataStructTemp = dataStruct;
            dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
            mnsSubDataFolder = [mnsDataFolder 'perm_' num2str(i) '_nl_' num2str(nl) '\'];
            [breakingReactionsPerm(:,i), noOfBreakingPermutations(i)] = ...
                permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsFolder, mnsSubDataFolder);
        end
    end
    mnsScanResults.p.pMode = 'only max no of reactions';
    mnsScanResults.p.breakingReactions = breakingReactionsPerm;
    mnsScanResults.p.noMaxBreakingPermutations = noOfBreakingPermutations;
    mnsScanResults.p.pvalueMaxBreaking = length(find(noOfBreakingPermutations > 0))/initMNS.permutations;
    if(initMNS.verboseScan == 2)
        disp('Permutations done')
    end
end

%% calc p-values for the complete scan ranges
if isfield(initMNS ,'determinePvalue') && initMNS.determinePvalue == 2
    if(initMNS.verboseScan == 2)
        disp('Start permutations: scan mode')
    end
    h = waitbar(0, ['0 of ' num2str(initMNS.permutations) ' permutations done'],'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0)
    waitObject = onCleanup(@() delete(h));
    initMnsTemp = initMNS;
    freqMaxBreakingPermutations = zeros(initMNS.permutations,length(nlRange)-1);
    freqSumBreakingReactions = zeros(initMNS.permutations,length(nlRange)-1);
    freqNoBreakingReactions = zeros(initMNS.permutations,length(nlRange)-1);
    freqMaxSumBreakingReactions = zeros(initMNS.permutations,(length(nlRange)-1)^2);
    modelTemp.iaMat = model.iaMat;
    modelTemp.metaboliteId = model.metaboliteId;
    dataStructTemp = dataStruct;
    
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

    copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'), fullfile(mnsDataFolder, 'nCliqueNoPerm.txt'));
    if ~isempty(mnsDataFolderCopy)
        copyfile(fullfile(mnsBaseFolder, 'opengm', mnsExec),...
            fullfile(mnsDataFolderCopy, mnsExec));
        copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'),...
            fullfile(mnsDataFolderCopy, 'nCliqueNoPerm.txt'));
        if(initMNS.input.edgeCliqueFileExists)
            copyfile(initMNS.input.edgeCliqueFile,...
            fullfile(mnsDataFolderCopy, 'cliques.txt'));
            initMnsTemp.input.edgeCliqueFile = fullfile(mnsDataFolderP, 'cliques.txt');
        end
    end
    if isfield(initMNS, 'parallel') && initMNS.parallel.useCluster
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
        elseif ~strcmp(poolobj.Cluster.Profile, initMNS.parallel.clusterName) || initMNS.parallel.cores ~= poolobj.NumWorkers
            poolobj.delete();
            poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
        end
    end
    for i = 1:initMnsTemp.permutations
        dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
        noOfBreakingPermutationsTemp = zeros(1,length(nlRange)-1);
        %%
        if isfield(initMNS, 'parallel') && initMNS.parallel.useCluster
            parfor j = 1:length(nlRange)-1
                mnsSubDataFolder = fullfile(mnsDataFolderP,[ 'perm_' num2str(i) '_nl_' num2str(nlRange(j))]);
                [breakingReactionsPerm(:,j), noOfBreakingPermutationsTemp(j)] = ...
                    mns_permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsDataFolderP, mnsSubDataFolder, nlRange(j),mnsExecFolderP);
            end
        else
            for j = 1:length(nlRange)-1
                mnsSubDataFolder = fullfile(mnsDataFolderP,[ 'perm_' num2str(i) '_nl_' num2str(nlRange(j))]);
                [breakingReactionsPerm(:,j), noOfBreakingPermutationsTemp(j)] = ...
                    mns_permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsFolder, mnsSubDataFolder, nlRange(j),mnsExecFolderP);
            end
        end
        %%
        sumBreaking = sum(breakingReactionsPerm,2);
        maxBreaking = zeros(size(sumBreaking));
        for j = 1:length(maxBreaking)
            idx = find(breakingReactionsPerm(j,:) == 1);
            if ~isempty(idx)
                maxBreaking(j) = max(idx(end));
            end
            
        end
        
        maxSumBreaking = sumBreaking.*maxBreaking;
        freqMaxBreakingPermutations(i,:) = hist(maxBreaking, 1:length(nlRange)-1);
        freqSumBreakingReactions(i,:) = hist(sumBreaking, 1:length(nlRange)-1);
        freqMaxSumBreakingReactions(i,:) = hist(maxSumBreaking, 1:(length(mnsScanResults.nlRange)-1)^2);
        freqNoBreakingReactions(i,:) = noOfBreakingPermutationsTemp;
        
        tp2 = toc(t1);
        dt = (tp2-tp1)/i;
        hours = floor(dt*(initMnsTemp.permutations-i)/3600);
        mins = round((dt*(initMnsTemp.permutations-i)-hours*3600)/60);
%         waitbar(i/initMNS.permutations, h, [{[num2str(i) ' of ' num2str(initMNS.permutations) ' permutations done']}; {['Time remaining: ' num2str(hours) ' h ' num2str(mins) ' min ']}]);
        waitbar(i/initMNS.permutations, h, ...
            [num2str(i) ' of ' num2str(initMNS.permutations) ' permutations done; '  num2str(hours) ' h ' num2str(mins) ' min left']);
        if getappdata(h,'canceling')
            delete(h);
            mnsScanResults = [];
            return
        end
    end
    delete(h);
%     poolobj.delete;
    
    mnsScanResults.p.pMode = 'scan mode';
    mnsScanResults.p.freqMaxBreakingReactions = freqMaxBreakingPermutations;
    mnsScanResults.p.pvalueMaxBreaking = zeros(1,length(nlRange)-1);
    mnsScanResults.p.freqSumBreakingReactions = freqSumBreakingReactions;
    mnsScanResults.p.pvalueSumBreaking = zeros(1,length(nlRange)-1);
    mnsScanResults.p.freqMaxSumBreakingReactions = freqMaxSumBreakingReactions;
    mnsScanResults.p.pvalueMaxSumBreaking = zeros(1,(length(mnsScanResults.nlRange)-1)^2);
    mnsScanResults.p.freqNoBreakingReactions = freqNoBreakingReactions;
    mnsScanResults.p.pvalueNoBreaking = zeros(1,length(nlRange)-1);   
    
    for i = length(mnsScanResults.nlRange)-1:-1:1
        mnsScanResults.p.pvalueMaxBreaking(i) = length(find(sum(freqMaxBreakingPermutations(:,i:end),2) > 0))/initMNS.permutations;
        mnsScanResults.p.pvalueSumBreaking(i) = length(find(sum(freqSumBreakingReactions(:,i:end),2) > 0))/initMNS.permutations;
        mnsScanResults.p.pvalueNoBreaking(i) = length(find(sum(freqNoBreakingReactions(:,i:end),2) > 0))/initMNS.permutations;
    end
    for i = (length(mnsScanResults.nlRange)-1)^2:-1:1
        mnsScanResults.p.pvalueMaxSumBreaking(i) = length(find(sum(mnsScanResults.p.freqMaxSumBreakingReactions(:,i:end),2) > 0))/initMNS.permutations;
    end
    if(initMNS.verboseScan == 2)
        disp('Permutations done')
    end
end

%% plot the neighborhood fracture and module labels over time
if plotResults
    figure, imagesc(mnsScanResults.breakingReaction)
    for i = 1:size(mns.results.reactionPairs,1)
        if isfield(model, 'mat')
            idxRp = find(model.mat(:,mns.results.reactionPairs(i,1)) ~= 0 & model.mat(:,mns.results.reactionPairs(i,2)) ~= 0);
            label{i} = [model.rpId{idxRp} ': ' model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}];
        else
            label{i} = [model.metaboliteName{mns.results.reactionPairs(i,1)} ' <-> ' model.metaboliteName{mns.results.reactionPairs(i,2)}];
        end
    end
    mnsScanResults.reactionLabels = label;
    colormap([1 1 1; 0.8 0 0.2])
    title('Neighborhood fractures')
    set(gca, 'YTick', 1:size(mns.results.reactionPairs,1), 'YTickLabel', label);
    set(gca, 'XTick', 1:length(nlRange), 'XTickLabel', nlRange);
    xlabel('\lambda_1')
    
    figure, imagesc(mnsScanResults.clustLabel)
    title('Module labels')
    set(gca, 'YTick', 1:size(mns.results.reactionPairs,1), 'YTickLabel', model.metaboliteName);
    set(gca, 'XTick', 1:length(nlRange), 'XTickLabel', nlRange);
    xlabel('\lambda_1')
    cmap = colormap(jet(mnsScanResults.mns{1}.parameters.noOfLabels));
    colormap(cmap)
    xmin = floor(min(dataStruct.data));
    xmax = ceil(max(dataStruct.data));
    xVec = xmin:0.1:xmax;
    
    % plot observation potential distributions
%     figure, hold on
%     for i = 1:mnsScanResults.mns{1}.parameters.noOfLabels
%         prob(:,i) = exp( - (xVec-mnsScanResults.mns{1}.results.groupMean(i)).^2 / (2*mnsScanResults.mns{1}.results.groupStd(i)^2));
%         plot(xVec, prob(:,i), '-', 'LineWidth', 2, 'Color', cmap(i,:));
%     end
    
    
    % plot last breaking l1 vs sum of breaking reactions
%     sumBreaking = sum(mnsScanResults.breakingReaction,2);
%     maxBreaking = zeros(size(sumBreaking));
%     for i = 1:length(maxBreaking)
%         idx = find(mnsScanResults.breakingReaction(i,:) == 1);
%         if ~isempty(idx)
%             maxBreaking(i) = max(idx(end));
%         end
%     end
%     a = rand(length(sumBreaking),1);
%     a(a < 0.5) = -0.25;
%     a(a >= 0.5) = 0.25;
%     
%     b = rand(length(sumBreaking),1);
%     b(b < 0.5) = -0.25;
%     b(b >= 0.5) = 0.25;
%     
% %     b = rand(length(sumBreaking),1);
%     figure,
%     h = plot(sumBreaking+a, maxBreaking+b, 'o');
%     xlabel('# Breaking reactions');
%     ylabel('Max. Breaking reactions');
%     customDataCursor(h, label)
    
    
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
mnsScanResults.combinedParameterizations = false;
end
%%

function absDiffSubProd = getAbsDiffSubProdReaction(reactionPairs,data)
absDiffSubProd = zeros(size(reactionPairs,1),1);

for i = 1:length(absDiffSubProd)
    absDiffSubProd(i) = abs(data(reactionPairs(i,1))-data(reactionPairs(i,2)));
end

end