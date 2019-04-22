%% MNS example for univariate data with single parameterizations

%% Summary:
% Example code to infer sites of metabolic regulation in fibroblasts with transketolase
% knockdowns (data from Kuehne et al, Mol Cell. 2015 Aug 6;59(3):359-71.)

%% load the data and model
load('mns_example_univariate_mns_fibroblasts_TK_knockdown - WS.mat')
% the data structure dataStrucTK contains log2(FC) data comparing
% fibroblasts with transketolase knockdown and f
%% initialize mns
%go to mns folder and exectute
mns_initialize

%% run MNS inference and calculate p-values
% 
data = dataStructTK;
model = KEGG_HSA_MNS_redGlycPPP;
model.metaboliteName = metAbbreviations;

% initalize MNS parameterization
% change number of hidden states to 5
% calculate the p-value
initMNS = mns_generateInitMNS('noOfclusters', 5, 'determinePvalue', 2); 

% run MNS inference 
mnsScanResultsTKKD = mns_scan2state(model,data, 'TKKD', 0, initMNS, true);

% change colormap of module label distribution
cmap = [236 0 140; 46 49 146; 0 174 239; 0 166 81;255 242 0]/255;
colormap(cmap)

clear cmap data initMNS model
%% generate outtable of the inference results sorted according rankSum
tkkd_results_table = mns_scanResult2table(mnsScanResultsTKKD, model, 'rankSum');

%% find the rank of a certain gene: TKT (Transketolase)
[rankSum, rankMax] = mns_scanFindGeneRank(mnsScanResultsTKKD, model,'TKT', 1);
disp(['Best rank of TKT according to total number of fractures: ' num2str(rankSum(1))]);

%% Show the metabolite changes arround the first ranked reactions in a biograph
mns_scanResults2biograph(model, mnsScanResultsTKKD, 'rankSum', 1, 6)

%% save the results in excel worksheets for import to cytoscape
% http://www.cytoscape.org/

model = KEGG_HSA_MNS_redGlycPPP;
model.metaboliteName = metAbbreviations;
mns_scan2state2cytoscape(model, mnsScanResultsTKKD, 'MNS Example - TKKD')

clear model