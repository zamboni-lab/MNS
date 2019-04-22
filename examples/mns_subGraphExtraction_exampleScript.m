% example: how to use the metabolic network segmentation site inference
% tool for sub graph extraction on comparative metabolomics data
%% initialize the MNS toolbox
mns_initialize

%% load the model in this case KEGG_HSA_MNS
load KEGG_HSA_MNS
model = KEGG_HSA_MNS;

%% load the data
load mns_subGraphExtraction_exampleWS
%% generate data structure
% required fiaExp and diff structure

nameTag = 'H2O2vsWT';
% select the data type of your annotation
% fiaMiner v2.0 -> 'fiaExp'
% fiaMiner v3.0 -> 'fiaExp v3.0'
% One column vector with metabolite id for each measurement -> 'MetIdList'
mnsDataStruct.dataType = 'fiaExp';
% add significance data (one column) to your structure. Data could be
% -log10(pval), abs(z-scores) etc...
mnsDataStruct.data = -log10(data.pval);
mnsDataStruct.annotation = data.annotation;
%% initialize all the parameters for the scanning

% since subgraph extraction mean values depend heavily on the min and max
% data this information is needed
temp.minData = min(mnsDataStruct.data);
temp.maxData = max(mnsDataStruct.data);
% if you want to change the std scanning range you can change/add the following
% temp.std0 = [2 3 4 5];
% temp.std1 = [-1 0 1 2];
% this creates a structure with the default parameter sets needed for the
% subgraph extraction
initMNS = mns_generateInitMNS('modeSubGraphExtraction', temp);
% if you want to change a parameter you can do it like that
% initMNS = mns_generateInitMNS('nL1steps',20);
initMNS.subGraphExtraction.mergeFreqCo = 0.2
%% start scanning
scanRange = 0; %defines range of the nL1 scanning, if 0 nL1 gets determined automatically
mnsSubGraphsResults = mns_extractSubGraphs(model,mnsDataStruct, nameTag, scanRange,initMNS);

%% subGraphs 2 outtable
% Rank is defined here as total number of breaks

outTable = mns_subGraphResult2table(mnsSubGraphsResults, model);
%% show sub graph results as biograph
subGraphIdx = [2]% 9 11 13 14];
for i = 1:length(subGraphIdx)
mns_subGraphResult2bioGraph(mnsSubGraphsResults.mergedSubGraphStruct, mnsDataStruct.data,...
     data.fc_log2, mnsSubGraphsResults.mns{1},model.metaboliteName, subGraphIdx(i))
end
%% mnsResults2cytoscape
% generate two xls-files with edge and node information that can imported
% into cytoscape. In the mns\external\cytoscapeStyle\ path you can find
% styles to visualize and explore the results.
subGraphIdx = 1;
mns_subGraphResult2cytoscape(mnsSubGraphsResults.mergedSubGraphStruct, mnsDataStruct.data,...
     data.fc_log2, mnsSubGraphsResults.mns{1},model.metaboliteName,nameTag, subGraphIdx)
