function [edgeStruct, edgeHead, nodeStruct, nodeHeader] = mns_subGraphResult2cytoscape(subGraphStruct, pvalData,fcData,mns,metName,nameTag,idx,save2File)
if nargin < 6 || isempty(nameTag)
    nameTag = 'tempSubGraph';
end
if nargin < 7
    idx = 1;
end
if nargin < 8
    save2File = true;
end
%% make edge file
metIdx = mns.observations.metId;
edges = subGraphStruct.subGraphDetails(idx).edges;
edgeFrequency = subGraphStruct.subGraphDetails(idx).edgeIdxFrequency;
noOfEdges = size(subGraphStruct.subGraphDetails(idx).edges,1);
edgeHead = {'SubId','EdgeId','ProdId','Frequency'};
edgeStruct = cell(noOfEdges,length(edgeHead));
for i = 1:noOfEdges
    edgeStruct(i,:) = [metIdx(edges(i,1)) {[metIdx{edges(i,2)} '-' metIdx{edges(i,2)}]}...
        metIdx(edges(i,2)) num2cell(edgeFrequency(i))];
end
if save2File
    xlswrite([nameTag '_edges_' num2str(idx) '.xls'], [edgeHead;edgeStruct]);
end
%% make node file
nodeHeader = {'metId', 'metName', 'detected' 'log2fc', '-log10pval'};
uMetIdx = subGraphStruct.subGraphDetails(idx).metIdx;
noOfMets = length(uMetIdx);
nodeStruct = cell(noOfMets,length(nodeHeader));
model2dataIdx = mns.observations.iaId2dataId(uMetIdx);
for i = 1:noOfMets
    id = uMetIdx(i);
    if model2dataIdx(i) ~= 0
        nodeStruct(i,:) = [mns.observations.metId(id) metName(id) ...
            1 num2cell(fcData(model2dataIdx(i))) num2cell(-log10(pvalData(model2dataIdx(i))))];
    else
        nodeStruct(i,:) = [mns.observations.metId(id) metName(id) 0 0 0];
    end
end
if save2File
    xlswrite([nameTag '_nodes_' num2str(idx) '.xls'], [nodeHeader;nodeStruct]);
end
end