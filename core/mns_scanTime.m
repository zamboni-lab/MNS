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
function mnsScanResults = mns_scanTime(model,dataStruct, nameTag, nlRange, tlRange, initMNS, plotResults)
warning off

t1 = tic;
if nargin < 3 || isempty(nameTag)
    nameTag = 'temp';
end
if nargin < 4 || isempty(nlRange)
    nlRange = 0;
end
if nargin < 5 || isempty(nlRange)
    tlRange = 0;
end
if nargin < 6 || isempty(initMNS)
    initMNS = mns_generateInitMNS;
end
if nargin < 7 || isempty(plotResults)
    plotResults = false;
end
if isfield(dataStruct, 'dataReplicates')
    initMNS.biolReplicates = dataStruct.dataReplicates;
end
%% preparations
% remove all white spaces from the nametag
nameTag = regexprep(nameTag, '\s', '\_');

% get mns folder and create data Folder
mnsFolder = regexprep(which('mns_runTimeFrame.m'),'mns_runTimeFrame\.m','');
mnsBaseFolder = regexprep(which('mns_initialize.m'),'mns_initialize\.m','');
switch computer
    case 'PCWIN32'
        mnsExec = 'mns_timeFrame.exe';
    case 'PCWIN64'
        mnsExec = 'mns_timeFrame_64.exe';
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
    otherwise
        mnsDataFolder = [mnsDataFolderStem '/' nameTag '/'];
end

mkdir(mnsDataFolder)

%% coarse grained search for nl
%if nlRange = 0 , do coarse grained search for nl with tl = 0;
if nlRange == 0
    if(initMNS.verboseScan == 2)
        disp('Start coarse grained search for the range of nl...');
    end
    nlTemp = -10;
    initMNStemp = initMNS;
    initMNStemp.tL3min = 0;
    initMNStemp.tL3max = 0;
    
    i = 0;
    while(nlTemp <= 512)
        i = i+1;
        tempStem = [mnsDataFolder 'coarse_nl_' regexprep(sprintf('%.0E',2^nlTemp),'\-','m') '_tl_0'];
        mkdir(tempStem);
        copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(tempStem, 'nCliqueNoPerm.txt'));
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, fullfile(tempStem, 'cliques.txt'));
        end
        if(initMNS.verboseScan == 2)
            disp(2^nlTemp);
        end
        initMNStemp.nL1max = 2^nlTemp;
        initMNStemp.nL1min= 2^nlTemp;
        [mns] = mns_runTimeFrame(model, dataStruct, initMNStemp, tempStem, mnsExecFolder);
        if isempty(mns.results)
            mnsScanResults = [];
            return;
        end
        if i == 1
            initMNStemp.input.initLabels = mns.results.initLabels;
            initMNStemp.labeling = 'fix';
            initMNStemp.input.reactionPairs = mns.results.reactionPairs;
            initMNStemp.input.edgeCliqueFile = fullfile(tempStem, 'cliques.txt');
            initMNStemp.input.edgeCliqueFileExists = true;
            initMNStemp.input.nCliques = mns.nCliques;
            initMNS.input.edgeCliqueFile = fullfile(tempStem, 'cliques.txt');
            initMNS.input.edgeCliqueFileExists = true;
            initMNS.input.nCliques = mns.nCliques;
        end
        nlTemp = nlTemp+1;
        if nlTemp >= 10
%             nlRange = 0:2^nlTemp/initMNS.nL1steps:2^nlTemp;
            nlRange = [0 2.^(-10:(nlTemp+10)/(initMNS.nL1steps-1):nlTemp)];
%             nlRange = [0 2.^(-10:nlTemp)];
            if(initMNS.verboseScan == 2)
                disp('Finished coarsed grain search (nl)');
%                 disp(['nlRange: 0:' num2str(2^nlTemp/initMNS.nL1steps) ':' num2str(2^nlTemp)]);
            end
            break;
        end
        
        if sum(sum(mns.results.breakingReactions,1)) == 0
%             nlRange = 0:2^nlTemp/initMNS.nL1steps:2^nlTemp;
            nlRange = [0 2.^(-10:(nlTemp+10)/(initMNS.nL1steps-1):nlTemp)];
            if(initMNS.verboseScan == 2)
                disp('Finished coarsed grain search (nl)');
%                 disp(['nlRange: 0:' num2str(2^nlTemp/initMNS.nL1steps) ':' num2str(2^nlTemp)]);
            end
            break;
        end
        
    end
    clear mnsScanResults initMNStemp i mns;
    if(initMNS.verboseScan == 2)
        disp('Coarse grained search finished (nl)')
    end
end

%% coarse grained search for tl
% if tlRange = 0 , do coarse grained search for tl with nl = 0;
if tlRange == 0
    if(initMNS.verboseScan == 2)
        disp('Start coarse grained search for the range of tl...');
    end
    tlTemp = -10;
    initMNStemp = initMNS;
    initMNStemp.nL1min = 0;
    initMNStemp.nL1max = 0;
    
    i = 0;
    while(tlTemp <= 512)
        i = i+1;
        tempStem = [mnsDataFolder 'coarse_nl_0_tl_' num2str(tlTemp)];
        mkdir(tempStem);
        copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(tempStem, 'nCliqueNoPerm.txt'));
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, fullfile(tempStem, 'cliques.txt'));
        end
        if(initMNS.verboseScan == 2)
            disp(2^tlTemp);
        end
        initMNStemp.tL3max = 2^tlTemp;
        initMNStemp.tL3min= 2^tlTemp;
        [mns] = mns_runTimeFrame(model, dataStruct, initMNStemp, tempStem, mnsExecFolder);
        if isempty(mns.results)
            mnsScanResults = [];
            return;
        end
        if i == 1
            initMNStemp.input.initLabels = mns.results.initLabels;
            initMNStemp.labeling = 'fix';
            initMNStemp.input.reactionPairs = mns.results.reactionPairs;
            initMNStemp.input.edgeCliqueFile = fullfile(tempStem, 'cliques.txt');
            initMNStemp.input.edgeCliqueFileExists = true;
            initMNStemp.input.nCliques = mns.nCliques;
            initMNS.input.edgeCliqueFile = fullfile(tempStem, 'cliques.txt');
            initMNS.input.edgeCliqueFileExists = true;
            initMNS.input.nCliques = mns.nCliques;
        end
        tlTemp = tlTemp+1;
        if tlTemp >= 9
%             tlRange = 0:2^tlTemp/initMNS.tL3steps:2^tlTemp;
            tlRange = [0 2.^(-10:(tlTemp+10)/(initMNS.tL3steps-1):tlTemp)];
%             tlRange = [0 2.^(-10:tlTemp)];
            if(initMNS.verboseScan == 2)
%                 disp('Finished coarsed grain search (tl)');
%                 disp(['tlRange: 0:' num2str(2^tlTemp/initMNS.tL3steps) ':' num2str(2^tlTemp)]);
            end
            break;
        end
        
        if sum(sum(mns.results.breakingReactionTime,1)) == 0
%             tlRange = 0:2^tlTemp/initMNS.tL3steps:2^tlTemp;
            tlRange = [0 2.^(-10:(tlTemp+10)/(initMNS.tL3steps-1):tlTemp)];
%             tlRange = [0 2.^(-10:tlTemp)];
            if(initMNS.verboseScan == 2)
%                 disp('Finished coarsed grain search (tl)');
%                 disp(['tlRange: 0:' num2str(2^tlTemp/initMNS.tL3steps) ':' num2str(2^tlTemp)]);
            end
            break;
        end
        
    end
    clear mnsScanResults initMNStemp i mns;
    if(initMNS.verboseScan == 2)
        disp('Coarse grained search finished (tl)')
    end
end

%% start the scanning
if(initMNS.verboseScan == 2)
    disp('Start scanning...')
end
c = 0;
nTp = size(dataStruct.data,2);
mnsScanResults.amountBreakingReaction = zeros(length(nlRange), length(tlRange));
mnsScanResults.amountTemporalBreaks = zeros(length(nlRange), length(tlRange));
% scan through nlRange and tlRange sequentially
for i = 1:length(nlRange)
    flag = 0;
    for j = 1:length(tlRange)
        c = c+1;
        paraComb(c,:) = [i j];
        %make folder stem
%         tempStemFolder = [mnsDataFolder 'nl_' num2str(nlRange(i)) 'tl_' num2str(tlRange(j))];
        tempStemFolder = [mnsDataFolder 'nl_' num2str(i) '_tl_' num2str(j)];
        mkdir(tempStemFolder);
        copyfile(fullfile(clFolder, 'nCliqueNoPerm.txt'), fullfile(tempStemFolder, 'nCliqueNoPerm.txt'));
        if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
            copyfile(initMNS.input.edgeCliqueFile, fullfile(tempStemFolder, 'cliques.txt'));
        end
        if(initMNS.verboseScan == 2)
            disp('');
            disp(['Starting step ' num2str(c) ', nl = ' num2str(nlRange(i)) ', tl = ' num2str(tlRange(j))]);
            disp('');
        end
        initMNS.nL1max = nlRange(i);
        initMNS.nL1min= nlRange(i);
        initMNS.tL3max = tlRange(j);
        initMNS.tL3min= tlRange(j);
        [mns] = mns_runTimeFrame(model, dataStruct, initMNS, tempStemFolder, mnsExecFolder);
        if isempty(mns.results)
            mnsScanResults = [];
            return;
        end
        mnsTemp = mns;
        rmfield(mnsTemp, 'nCliques');
        if c == 1
            mnsScanResults.mns(c) = {mnsTemp};
        end
        mnsScanResults.breakingReaction(:,:,c) = mns.results.breakingReactions;
        tempidx = mns.observations.modelVarIdTimePoint+1;
%         tempidx = mnsScanResults.mns{c}.observations.modelVarIdTimePoint+1;
        mnsScanResults.clustLabel(:,:,c) = mns.results.inferedLabels(tempidx);
        mnsScanResults.temporalBreaks(:,:,c) = mns.results.breakingReactionTime;
        mnsScanResults.nL1(c) = nlRange(i);
        mnsScanResults.tL3(c) = tlRange(j);
        mnsScanResults.amountBreakingReaction(i,j) = sum(sum(mnsScanResults.breakingReaction(:,:,c)));
        mnsScanResults.amountTemporalBreaks(i,j) = sum(sum(mnsScanResults.temporalBreaks(:,:,c)));
        if i == 1 && j == 1
            initMNS.input.edgeCliqueFile = fullfile(tempStemFolder, 'cliques.txt');
            initMNS.input.edgeCliqueFileExists = true;
            initMNS.input.nCliques = mns.nCliques;
            initMNS.input.initLabels = mnsScanResults.mns{1,1}.results.initLabels;
            initMNS.labeling = 'fix';
            initMNS.input.reactionPairs = mnsScanResults.mns{1,1}.results.reactionPairs;
            mnsScanResults.maxNeighborBreakingReaction = zeros(size(mnsScanResults.breakingReaction,1),length(tlRange));
            mnsScanResults.sumNeighborBreakingReaction = zeros(size(mnsScanResults.breakingReaction,1),length(tlRange));
            mnsScanResults.maxTimeBreakingReaction = zeros(size(mnsScanResults.breakingReaction,1),length(nlRange));
            mnsScanResults.sumTimeBreakingReaction = zeros(size(mnsScanResults.breakingReaction,1),length(nlRange));
        end
        % save for maxReaction dependent on nl or tl state
        %% stil have to think about that
%         mnsScanResults.maxNeighborBreakingReaction(mnsScanResults.breakingReaction(:,:,c) == 1,j) = i;
%         mnsScanResults.maxTimeBreakingReaction(mnsScanResults.breakingReaction(:,i,j) == 1,i) = j;
    %%
        % if no moreBreakinReactions -> increase i or if j == 0 stop
        % scanning
        if mnsScanResults.amountBreakingReaction(i,j) == 0 && j == 1
%             if j == 1
%                 nlRange(i+1:end) = [];
%                 mnsScanResults.maxTimeBreakingReaction(:,i+1:end) = [];
%                 mnsScanResults.sumTimeBreakingReaction(:,i+1:end) = [];
%                 maxTimeId = max(max(mnsScanResults.maxTimeBreakingReaction))+1;
%                 tlRange(maxTimeId+1:end) = [];
%                 mnsScanResults.maxNeighborBreakingReaction(:,maxTimeId+1:end) = [];
%                 mnsScanResults.sumNeighborBreakingReaction(:,maxTimeId+1:end) = [];
%                 break;
%             else
%                 break;
%             end
            flag = 1;
            break;
        end
%         if sum(sum(mnsScanResults.breakingReaction(:,:,c),1)) == 0
%             break;
%         end
        if mnsScanResults.amountTemporalBreaks(i,j) == 0
            break;
        end
    end
    if flag
        break;
    end
%     if mnsScanResults.amountBreakingReaction(i,j) == 0
%         break;
%     end
end
%% find max of ranges
tlRangeMax = min(find(sum(mnsScanResults.amountTemporalBreaks,1) == 0));
nlRangeMax = min(find(sum(mnsScanResults.amountBreakingReaction,2) == 0));
if isempty(tlRangeMax)
    tlRangeMax = length(tlRange);
end
if isempty(nlRangeMax)
    nlRangeMax = length(nlRange);
end
mnsScanResults.nL1Range = nlRange(1:nlRangeMax);
mnsScanResults.tL3Range = tlRange(1:tlRangeMax);
mnsScanResults.amountTemporalBreaks = mnsScanResults.amountTemporalBreaks(1:nlRangeMax, 1:tlRangeMax);
mnsScanResults.amountBreakingReaction = mnsScanResults.amountBreakingReaction(1:nlRangeMax, 1:tlRangeMax);
% mnsScanResults.amountTemporalBreaks = zeros(length(nlRange), length(tlRange));
if(initMNS.verboseScan == 2)
    disp('Scanning finished')
end
% mnsScanResults.sumBreakingReaction = sum(mnsScanResults.breakingReaction,2);
% uMaxBreakingReaction = unique(mnsScanResults.maxBreakingReaction);
% 
% mnsScanResults.nlRange = nlRange;
% mnsScanResults.tlRange = tlRange;

% % calc p-value only for maxRecations
% if isfield(initMNS ,'determinePvalue') && initMNS.determinePvalue == 1
%     if(initMNS.verboseScan == 2)
%         disp('Start permutations: only max')
%     end
%     nl = nlRange(end-1);
% %     h = waitbar(0, ['0 of ' num2str(initMNS.permutations) ' permutations done']);
% %     dataStructTemp = dataStruct;
%     initMnsTemp = initMNS;
%     initMnsTemp.nL1max = nlRange(end-1);
%     initMnsTemp.nL1min = nlRange(end-1);
%     noOfBreakingPermutations = zeros(initMNS.permutations,1);
%     modelTemp.iaMat = model.iaMat;
%     modelTemp.metaboliteId = model.metaboliteId;
%     
%     if isfield(initMNS.parallel, 'useCluster') && initMNS.parallel.useCluster
%         poolobj = gcp('nocreate');
%         if isempty(poolobj) || ~strcmp(initMNS.parallel.clusterName, poolobj.Cluster.profile) || initMNS.parallel.cores ~= poolobj.Cluster.NumWorkers
%             poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
%         end
%         parfor i = 1:initMnsTemp.permutations
%             dataStructTemp = dataStruct;
%             dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
%             mnsSubDataFolder = [mnsDataFolder 'perm_' num2str(i) '_nl_' num2str(nl) '\'];
%             [breakingReactionsPerm(:,i), noOfBreakingPermutations(i)] = ...
%                 permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsFolder, mnsSubDataFolder);
%         end
%     else
%         for i = 1:initMnsTemp.permutations
%             dataStructTemp = dataStruct;
%             dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
%             mnsSubDataFolder = [mnsDataFolder 'perm_' num2str(i) '_nl_' num2str(nl) '\'];
%             [breakingReactionsPerm(:,i), noOfBreakingPermutations(i)] = ...
%                 permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsFolder, mnsSubDataFolder);
%         end
%     end
% %     delete(h);
% %     poolobj.delete;
%     mnsScanResults.p.pMode = 'only max no of reactions';
%     mnsScanResults.p.breakingReactions = breakingReactionsPerm;
%     mnsScanResults.p.noMaxBreakingPermutations = noOfBreakingPermutations;
%     mnsScanResults.p.pvalueMaxBreaking = length(find(noOfBreakingPermutations > 0))/initMNS.permutations;
%     if(initMNS.verboseScan == 2)
%         disp('Permutations done')
%     end
% end
% 
% % calc p-value only for maxRecations
% if isfield(initMNS ,'determinePvalue') && initMNS.determinePvalue == 2
%     if(initMNS.verboseScan == 2)
%         disp('Start permutations: scan mode')
%     end
%     h = waitbar(0, ['0 of ' num2str(initMNS.permutations) ' permutations done']);
%     initMnsTemp = initMNS;
%     freqMaxBreakingPermutations = zeros(initMNS.permutations,length(nlRange)-1);
%     freqSumBreakingReactions = zeros(initMNS.permutations,length(nlRange)-1);
%     freqNoBreakingReactions = zeros(initMNS.permutations,length(nlRange)-1);
%     diffMaxVsMeanBreakingReactionPermuted = zeros(initMNS.permutations,size(mnsScanResults.breakingReaction,1));
%     diffNextBreakingReactionPermuted = zeros(initMNS.permutations,size(mnsScanResults.breakingReaction,1));
%     modelTemp.iaMat = model.iaMat;
%     modelTemp.metaboliteId = model.metaboliteId;
%     dataStructTemp = dataStruct;
%     
%     tp1 = toc(t1);
%     if isfield(initMNS.parallel, 'parallelFolder') && ~isempty(initMNS.parallel.parallelFolder)
%         mnsDataFolderP = initMNS.parallel.parallelFolder;
%         mnsDataFolderCopy = initMNS.parallel.parallelFolderCopy;
%         mnsExecFolderP = mnsDataFolderP;
%     else
%         mnsDataFolderP = mnsDataFolder;
%         mnsDataFolderCopy = '';
%         mnsExecFolderP = '';
%     end
%     if ~isempty(mnsDataFolderCopy)
%         copyfile(fullfile(mnsFolder, 'opengm', mnsExec),...
%             fullfile(mnsDataFolderCopy, mnsExec));
%         copyfile(fullfile(mnsFolder, 'nCliqueNoPerm.txt'),...
%             fullfile(mnsDataFolderCopy, 'nCliqueNoPerm.txt'));
%         if(initMNS.input.edgeCliqueFileExists)
%             copyfile(initMNS.input.edgeCliqueFile,...
%             fullfile(mnsDataFolderCopy, 'cliques.txt'));
%             initMnsTemp.input.edgeCliqueFile = fullfile(mnsDataFolderP, 'cliques.txt');
%         end
%     end
%     for i = 1:initMnsTemp.permutations
%         dataStructTemp.data = dataStruct.data(randperm(size(dataStruct.data,1)));
%         noOfBreakingPermutationsTemp = zeros(1,length(nlRange)-1);
%         if isfield(initMNS, 'parallel') && initMNS.parallel.useCluster
%             poolobj = gcp('nocreate');
%             if isempty(poolobj)
%                 poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
%             elseif ~strcmp(poolobj.Cluster.Host, initMNS.parallel.clusterName) || initMNS.parallel.cores ~= poolobj.NumWorkers
%                 poolobj.delete();
%                 poolobj = parpool(initMNS.parallel.clusterName,initMNS.parallel.cores);
%             end
%             parfor j = 1:length(nlRange)-1
%                 mnsSubDataFolder = [mnsDataFolderP 'perm_' num2str(i) '_nl_' num2str(nlRange(j)) '\'];
%                 [breakingReactionsPerm(:,j), noOfBreakingPermutationsTemp(j)] = ...
%                     permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsDataFolderP, mnsSubDataFolder, nlRange(j),mnsExecFolderP);
%             end
%         else
%             for j = 1:length(nlRange)-1
%                 mnsSubDataFolder = [mnsDataFolderP 'perm_' num2str(i) '_nl_' num2str(nlRange(j)) '\'];
%                 [breakingReactionsPerm(:,j), noOfBreakingPermutationsTemp(j)] = ...
%                     permStartMNS(modelTemp, initMnsTemp, dataStructTemp, mnsFolder, mnsSubDataFolder, nlRange(j),mnsExecFolderP);
%             end
%         end
%         sumBreaking = sum(breakingReactionsPerm,2);
%         maxBreaking = zeros(size(sumBreaking));
%         for j = 1:length(maxBreaking)
%             idx = find(breakingReactionsPerm(j,:) == 1);
%             if ~isempty(idx)
%                 maxBreaking(j) = max(idx(end));
%             end
%             
%         end
%         diffMaxVsMeanBreakingReactionPermuted(i,:) = maxBreaking - mean(maxBreaking);
%         uMaxBreakingReaction = unique(maxBreaking);
%         for j = 2:length(uMaxBreakingReaction)
%             diffNextBreakingReactionPermuted(i,mnsScanResults.maxBreakingReaction == uMaxBreakingReaction(j)) ...
%                 = uMaxBreakingReaction(j) - uMaxBreakingReaction(j-1);
%         end
%         
%         freqMaxBreakingPermutations(i,:) = hist(maxBreaking, 1:length(nlRange)-1);
%         freqSumBreakingReactions(i,:) = hist(sumBreaking, 1:length(nlRange)-1);
%         freqNoBreakingReactions(i,:) = noOfBreakingPermutationsTemp;
%         
%         tp2 = toc(t1);
%         dt = (tp2-tp1)/i;
%         hours = floor(dt*(initMnsTemp.permutations-i)/3600);
%         mins = round((dt*(initMnsTemp.permutations-i)-hours*3600)/60);
%         waitbar(i/initMNS.permutations, h, [{[num2str(i) ' of ' num2str(initMNS.permutations) ' permutations done']}; {['Time remaining: ' num2str(hours) ' h ' num2str(mins) ' min ']}]);
%     end
%     delete(h);
% %     poolobj.delete;
%     
%     mnsScanResults.p.pMode = 'scan mode';
% %     mnsScanResults.p.breakingReactions = breakingReactionsPerm;
%     mnsScanResults.p.freqMaxBreakingReactions = freqMaxBreakingPermutations;
%     mnsScanResults.p.pvalueMaxBreaking = zeros(1,length(nlRange)-1);
%     mnsScanResults.p.freqSumBreakingReactions = freqSumBreakingReactions;
%     mnsScanResults.p.pvalueSumBreaking = zeros(1,length(nlRange)-1);
%     mnsScanResults.p.freqNoBreakingReactions = freqNoBreakingReactions;
%     mnsScanResults.p.pvalueNoBreaking = zeros(1,length(nlRange)-1);
%     mnsScanResults.p.pvalueDiffMaxVsMeanBreaking = zeros(size(mnsScanResults.breakingReaction,1),1);
%     mnsScanResults.p.pvalueDiffMaxVsNextBreaking = zeros(size(mnsScanResults.breakingReaction,1),1);
%     
%     
%     for i = length(mnsScanResults.nlRange)-1:-1:1
%         mnsScanResults.p.pvalueMaxBreaking(i) = length(find(sum(freqMaxBreakingPermutations(:,i:end),2) > 0))/initMNS.permutations;
% %         mnsScanResults.p.pvalueMaxBreaking(i) = sum(sum(freqMaxBreakingPermutations(:,i:end),2))/(size(mnsScanResults.breakingReaction,1)*initMNS.permutations);
%         mnsScanResults.p.pvalueSumBreaking(i) = length(find(sum(freqSumBreakingReactions(:,i:end),2) > 0))/initMNS.permutations;
% %         mnsScanResults.p.pvalueSumBreaking(i) = sum(sum(freqSumBreakingReactions(:,i:end),2))/(size(mnsScanResults.breakingReaction,1)*initMNS.permutations);
%         mnsScanResults.p.pvalueNoBreaking(i) = length(find(sum(freqNoBreakingReactions(:,i:end),2) > 0))/initMNS.permutations;
%     end
%     for i = 1:length(mnsScanResults.maxBreakingReaction)
%         mnsScanResults.p.pvalueDiffMaxVsMeanBreaking(i) = ...
%             length(find(diffMaxVsMeanBreakingReactionPermuted >= mnsScanResults.diffMaxVsMeanBreakingReaction(i)))/(size(diffMaxVsMeanBreakingReactionPermuted,1)*size(diffMaxVsMeanBreakingReactionPermuted,2));
%         mnsScanResults.p.pvalueDiffMaxVsNextBreaking(i) = ...
%             length(find(diffNextBreakingReactionPermuted >= mnsScanResults.diffNextBreakingReaction(i)))/(size(diffNextBreakingReactionPermuted,1)*size(diffNextBreakingReactionPermuted,2));
%     end
%     if(initMNS.verboseScan == 2)
%         disp('Permutations done')
%     end
% end


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
    colormap([1 1 1; 0.75 0 0])
    title('Breaking reactions')
    set(gca, 'YTick', 1:size(mns.results.reactionPairs,1), 'YTickLabel', label);
    set(gca, 'XTick', 1:length(nlRange), 'XTickLabel', nlRange);
    xlabel('\lambda_1')
    
    figure, imagesc(mnsScanResults.clustLabel)
    title('ClusterLabels')
    set(gca, 'YTick', 1:size(mns.results.reactionPairs,1), 'YTickLabel', model.metaboliteName);
    set(gca, 'XTick', 1:length(nlRange), 'XTickLabel', nlRange);
    xlabel('\lambda_1')
    cmap = colormap(jet(mnsScanResults.mns{1}.parameters.noOfLabels));
    cmap = [236 0 140; 46 49 146; 0 174 239; 0 166 81;255 242 0]/255;
    colormap(cmap)
    xmin = floor(min(dataStruct.data));
    xmax = ceil(max(dataStruct.data));
    xVec = xmin:0.1:xmax;
    xVec = -6:0.1:6;
    
    % plot observation potential distributions
    figure, hold on
    for i = 1:mnsScanResults.mns{1}.parameters.noOfLabels
        prob(:,i) = exp( - (xVec-mnsScanResults.mns{1}.results.groupMean(i)).^2 / (2*mnsScanResults.mns{1}.results.groupStd(i)^2));
        plot(xVec, prob(:,i), '-', 'LineWidth', 2, 'Color', cmap(i,:));
    end
    
    
    % plot last breaking l1 vs sum of breaking reactions
    sumBreaking = sum(mnsScanResults.breakingReaction,2);
    maxBreaking = zeros(size(sumBreaking));
    for i = 1:length(maxBreaking)
        idx = find(mnsScanResults.breakingReaction(i,:) == 1);
        if ~isempty(idx)
            maxBreaking(i) = max(idx(end));
        end
    end
    a = rand(length(sumBreaking),1);
    a(a < 0.5) = -0.25;
    a(a >= 0.5) = 0.25;
    
    b = rand(length(sumBreaking),1);
    b(b < 0.5) = -0.25;
    b(b >= 0.5) = 0.25;
    
    b = rand(length(sumBreaking),1);
    figure,
    h = plot(sumBreaking+a, maxBreaking+b, 'o');
    xlabel('# Breaking reactions');
    ylabel('Max. Breaking reactions');
    customDataCursor(h, label)
    
    
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

function [breakingReactions, noOfBreakingPermutations] = permStartMNS(model, initMNS, dataStruct, mnsFolder, mnsSubDataFolder, nl,mnsExecFolderP)
    if nargin > 5 && ~isempty(nl)
        initMNS.nL1max = nl;
        initMNS.nL1min = nl;
    end
    
    mkdir(mnsSubDataFolder);

    copyfile([mnsFolder 'nCliqueNoPerm.txt'], [mnsSubDataFolder 'nCliqueNoPerm.txt']);
    if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
        copyfile(initMNS.input.edgeCliqueFile, [mnsSubDataFolder 'cliques.txt']);
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

function [mnsCell,breakingReactions, noOfBreakingPermutations] = startMNSparallel(model, initMNS, dataStruct, mnsFolder, mnsSubDataFolder, nl)
    if nargin > 5 && ~isempty(nl)
        initMNS.nL1max = nl;
        initMNS.nL1min = nl;
    end
    
    mkdir(mnsSubDataFolder);

    copyfile([mnsFolder 'nCliqueNoPerm.txt'], [mnsSubDataFolder 'nCliqueNoPerm.txt']);
    if isfield(initMNS.input, 'edgeCliqueFileExists') && initMNS.input.edgeCliqueFileExists
        copyfile(initMNS.input.edgeCliqueFile, [mnsSubDataFolder 'cliques.txt']);
    end
    [mnsTEMP, ~, ~] = mns_run2state(model, dataStruct, initMNS, mnsSubDataFolder);
    breakingReactions = mnsTEMP.results.breakingReactions;
    noOfBreakingPermutations = sum(mnsTEMP.results.breakingReactions);
    mnsCell = {mnsTEMP};
    try
        rmdir(mnsSubDataFolder,'s');
    catch err
        disp(err)
    end
end
