% example: how to use the metabolic network segmentation / site inference
% tool for two state data
%% initialize the MNS toolbox
mns_initialize

%% load the model in this case KEGG_ECO_MNS
load KEGG_ECO_MNS
model = KEGG_ECO_MNS;

%% load the data
load mns_mnsScan2state_exampleWS
%% generate data structure

nameTag = 'purMvsWT';
% select the data type of your annotation
% fiaMiner v2.0 -> 'fiaExp'
% fiaMiner v3.0 -> 'fiaExp v3.0'
% One column vector with metabolite id for each measurement -> 'MetIdList'
mnsDataStruct.dataType = 'fiaExp';
% add comparative data (one column) to your structure. Data could ne
% log2(fc), z-scores etc...
mnsDataStruct.data = data.fc_log2;
mnsDataStruct.annotation = data.annotation;
%% initialize all the parameters for the scanning
% this creates a structure with the default parameter sets
initMNS = mns_generateInitMNS%('meanType', 'initLabels', 'noOfclusters', 5);

% if you want to change a parameter you can do it like that
% initMNS = mns_generateInitMNS('nL1steps',20);

%% start scanning
scanRange = 0; %defines range of the nL1 scanning, if 0 nL1 gets determined automatically
mnsScanResults = mns_scan2state(model,mnsDataStruct, nameTag, scanRange, initMNS);

%% mnsResults 2 outtable
% Rank is defined here as total number of breaks
sortBy = 'rankMax'; %Sorting of the outable; options rankMax, rankSum, rankMaxRankSum
outTable = mns_scanResult2table(mnsScanResults, model,sortBy);

%% mnsResults2cytoscape
% generate two xls-files with edge and node information that can imported
% into cytoscape. In the mns\external\cytoscapeStyle\ path you can find
% styles to visualize and explore the results.
mns_scan2state2cytoscape(model, mnsScanResults, nameTag)

%% find gene rank
% k = amount of neighbouring considered for ranknumber of neighbour of gene
pert = 'purM';
k = 2;
[rankSum, rankMax] = mns_scanFindGeneRank(mnsScanResults, model,pert, k);

% Description of output
% col 1 is always for gene(reaction) itself, 2nd column is for the best
% rank of the 1st closest reactions, 3rd columns for the best rank of the
% 2nd closest reactions etc
% output arrays
% rankSum = tiedrank of the reaction for the total number of breaks
% rankMax = tiedrank of the reaction for the max breaks/fractures
% rankUSum = unique rank of the reaction for the total number of breaks
% rankUMax = unique rank of the reaction for the max breaks
% maxRankUSum = what is the total number of different unique ranks for
% total number of breaks
% maxRankUMax = what is the total number of different unique ranks for
% max breaks
% Difference tiedrank unique rank: Example
% TotalBreaks    tiedrank     uniqueRank
%     10          -> 1.5        -> 1
%     10          -> 1.5        -> 1
%     8           -> 5          -> 2
%     8           -> 5          -> 2
%     8           -> 5          -> 2
%     8           -> 5          -> 2
%     8           -> 5          -> 2
%     6           -> 8          -> 3
