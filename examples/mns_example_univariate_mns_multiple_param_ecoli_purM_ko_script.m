%% Example 2
% Identification of regulatory sites in E. coli with purM knockout using
% MNS for univariate data with a single parameterization

%% Summary:
% This code examplifies how to use the MNS inference for sites of metabolic
% regulations with multiple model paramterizations

%% initialize MNS
% go to mns folder and perform the execute command:
mns_initialize

%% load the data
% the datastructure contains log2(FC) data comparing E. coli + glucose with
% purM KO and wild type E. coli both cultured in M9 minimal medium
load('mns_example_univariate_mns_multiple_param_ecoli_purM_ko - WS')

%% load the model
% KEGG Ecoli
load KEGG_ECO_MNS
model = KEGG_ECO_MNS;

%% initialize different MNS parameterizations

% generate initMNS_3
% (parameter combination 3)
initMNS_2 = mns_generateInitMNS('verboseScan', 0, 'stdType', 'one group');

% generate initMNS_3
% (parameter combination 3)
initMNS_3 = mns_generateInitMNS('verboseScan', 0, 'stdType', 'one group - all data - factor'...
    ,'stdVal',1);

%% run MNS inference for multiple parameterizations
mnsScanResults_p2_3 = mns_scan2state_multipleParameterizations(model,dataStruct_purM_KO,'purM_KO', 0, initMNS_2, initMNS_3);

%% find gene rank of purM KO
[rankSum, rankMax] = mns_scanFindGeneRank(mnsScanResults_p2_3,model, 'purM', 1);
disp(['Best rank of purM KO according to rankproduct of total number of fractures: ' num2str(rankSum(1))]);

%% generate results table
resultsTable = mns_scanResult2table(mnsScanResults_p2_3, model,'rankproductSum');