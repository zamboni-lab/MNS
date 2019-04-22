function [edgeStruct, edgeHead, nodeStruct, nodeHeader] = mns_subGraphMetGeneResult2cytoscape(subGraphStruct, pvalData,fcData,...
    pvalGeneData,fcGeneData,mns,metName,geneName,rpName,nameTag, idx, collapseGeneIds,ws, save2file)
if nargin < 10 || isempty(nameTag)
    nameTag = 'tempSubGraph';
end
if nargin < 11 || isempty(idx)
    idx = 1;
end
if nargin < 12
    collapseGeneIds = false;
end
if nargin < 13
    ws = num2str(idx);
end
if nargin < 14
    save2file = true;
end
%% make edge file
metIdx = mns.observations.met.metId;
uCmpdIdx = subGraphStruct.subGraphDetails(idx).cmpdIdx;
noOfCmpds = length(uCmpdIdx);
edges = subGraphStruct.subGraphDetails(idx).edges;
edgeFrequency = subGraphStruct.subGraphDetails(idx).edgeIdxFrequency;
noOfEdges = size(subGraphStruct.subGraphDetails(idx).edges,1);
edgeHead = {'SubId','EdgeId','ProdId','Frequency'};
edgeStruct = cell(noOfEdges,length(edgeHead));
spontCount = 0;
for i = 1:noOfCmpds
    temp = mns.observations.iaId2dataId{uCmpdIdx(i)};
    if temp ~= 0
        
    tempName = unique(geneName(temp));
    tempStr = tempName{1};
    for j = 2:length(tempName)
        tempStr = [tempStr '; ' tempName{j}];
    end
    tempGeneName(uCmpdIdx(i)) = {tempStr};
    else
        spontCount = spontCount +1;
        tempGeneName(uCmpdIdx(i)) = {['spontanous_' num2str(spontCount)]};
    end
end
if collapseGeneIds
    tempName = tempGeneName;
else
    tempName = rpName;
end
for i = 1:noOfEdges
    if mns.observations.isMet(edges(i,1))
        subName = metIdx(edges(i,1)-mns.observations.noOfRp);
    elseif mns.observations.isRp(edges(i,1))
        subName = tempName(edges(i,1));
    end
    if mns.observations.isMet(edges(i,2))
        prodName = metIdx(edges(i,2)-mns.observations.noOfRp);
    elseif mns.observations.isRp(edges(i,2))
        prodName = tempName(edges(i,2));
    end
    edgeStruct(i,:) = [subName {[subName{:} '-' prodName{:}]}...
        prodName num2cell(edgeFrequency(i))];
end
[~,idxEdge] = unique(edgeStruct(:,2));
if save2file
    xlswrite([nameTag '_edges_' ws '.xls'], [edgeHead;edgeStruct(idxEdge,:)]);
end

%% make node file
nodeHeader = {'metId', 'metName', 'detected' 'log2fc', '-log10pval' 'type'};
% uMetIdx = subGraphStruct.subGraphDetails(idx).metIdx;

nodeStruct = cell(noOfCmpds,length(nodeHeader));
model2dataIdx = mns.observations.iaId2dataId(uCmpdIdx);
detected = mns.observations.detected(uCmpdIdx);
for i = 1:noOfCmpds
    id = uCmpdIdx(i);
    if mns.observations.isMet(id)
        if model2dataIdx{i} ~= 0
            nodeStruct(i,:) = [metIdx(id-mns.observations.noOfRp) metName(id-mns.observations.noOfRp) ...
                1 num2cell(fcData(model2dataIdx{i}))...
                num2cell(-log10(pvalData(model2dataIdx{i}))) ...
                'metabolite'];
        else
            nodeStruct(i,:) = [metIdx(id-mns.observations.noOfRp) metName(id-mns.observations.noOfRp) 0 0 0 'metabolite'];
        end
    elseif mns.observations.isRp(id)
        if detected(i) == 1
            pval = mns.observations.data(id);
            tempFc = fcGeneData(mns.observations.iaId2dataId{id});
            fc = median(tempFc);
            nodeStruct(i,:) = [tempName(id) tempGeneName(id) 1 fc pval 'gene'];
        else
            nodeStruct(i,:) = [tempName(id) tempGeneName(id) 0 0 0 'gene'];
        end
    end
end
[~,idxNode] = unique(nodeStruct(:,1));
if save2file
    xlswrite([nameTag '_nodes_' ws '.xls'], [nodeHeader;nodeStruct(idxNode,:)]);
end



end