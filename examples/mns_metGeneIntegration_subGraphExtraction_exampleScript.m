%example script for the subgraph extraction using metabolomics and
%transcriptomics data using the mns toolbox
%% initialize the MNS toolbox
mns_initialize
%% load model and generate bipartite model
% load the model
load KEGG_HSA_MNS
model = KEGG_HSA_MNS;

%generate for each gene one reaction pair
splitGeneRp = true;
model_bipartite = mns_rpair2bipartiteModel(model, splitGeneRp);
%% load the data
load mns_metGeneIntegration_subGraphExtraction_exampleWS

%% generate the data structure
nameTag = 'HiCa168';

dataStruct.dataType = 'fiaExp';
% select the data type of your annotation
% fiaMiner v2.0 -> 'fiaExp'
% fiaMiner v3.0 -> 'fiaExp v3.0'
% One column vector with metabolite id for each measurement -> 'MetIdList'

% significane data for metabolites
dataStruct.data = -log10(data.pval);
%annotation of metabolites
dataStruct.annotation = data.annotation;
% select gene annotation (so far only gene symbol is implemented)
dataStruct.geneDataType = 'GeneSymbol';
% significane data for genes
dataStruct.geneData = -log10(data.genePval);
% annotation for the genes
dataStruct.geneAnnotation = data.geneSymbol;

%% initialize the parameter set

% since subgraph extraction mean values depend heavily on the min and max
% data this information is needed
temp.minData = min(dataStruct.data);
temp.maxData = max(dataStruct.data);
temp.minGeneData = min(dataStruct.geneData);
temp.maxGeneData = max(dataStruct.geneData);
% if you want to change the std scanning range you can change/add the following
temp.geneStd0 = [9 12 15 18]./3;
temp.geneStd1 = [2 3 4 5];
temp.std0 = [2 2 4 4]./3;
temp.std1 = [0 0 1 1];
% this creates a structure with the default parameter sets needed for the
% subgraph extraction
initMNS = mns_generateInitMNS('modeSubGraphExtractionMetGene', temp);

% change the no of nL1 scanning steps to 20
initMNS.nL1steps = 20;

%change the inferenceParameter to 1 (default is two)
initMNS.input.inferenceParameter = 1;

%scan range; if 0 it is determined automatically
scanRange = 0;

mnsSubGraphsMetGene = mns_extractSubGraphs_MetGeneIntegration(model_bipartite,dataStruct, nameTag, scanRange,initMNS);
%% plot the subgraphs to biographs
subGraphIdx = 1:10;
for i = 1:length(subGraphIdx)
mns_subGraphMetGeneResult2bioGraph(mnsSubGraphsMetGene.mergedSubGraphStruct, data.pval,data.fc_log2,...
    data.genePval,data.geneFc_log2,mnsSubGraphsMetGene.mns{1},model_bipartite.metaboliteName,...
    dataStruct.geneAnnotation,model_bipartite.rpId,subGraphIdx(i),true)
end

%%  plot the subgraphs to cytoscape
subGraphIdx = [1:2 8];
for i = 1:length(subGraphIdx)
mns_subGraphMetGeneResult2cytoscape(mnsSubGraphsMetGene.mergedSubGraphStruct, data.pval,data.fc_log2,...
    data.genePval,data.geneFc_log2,mnsSubGraphsMetGene.mns{1},model_bipartite.metaboliteName,...
    dataStruct.geneAnnotation,model_bipartite.rpId,nameTag,subGraphIdx(i),true)
end

