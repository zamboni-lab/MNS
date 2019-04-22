function mns_subGraphMetGeneResult2bioGraph(subGraphStruct, pvalData,fcData,...
    pvalGeneData,fcGeneData,mns,metName,geneName,rpName,idx, collapseGeneIds,cmax)
if nargin < 10
    idx = 1;
end
if nargin < 11
    collapseGeneIds = true;
end
if nargin < 12
    cmax = 3;
end
%% make dualColorColormap

offset = 1;
upOffset = 1;

% cmax = 3;
cSize = cmax*10;
cmap = zeros(cSize,3);
cmap(offset+1:cSize,1) = 0:1/(cSize-offset-1):1;
if upOffset > 0
   cmap(cSize+1:cSize+offset,1) = 1; 
end

cmap2 = zeros(size(cmap,1)*2,3);
cmap2(1:size(cmap,1),3) = cmap(end:-1:1,1);
cmap2(size(cmap,1)+1:end,:) = cmap;
%% make dualColorColormap
N = 255;
X = [0.5: -1/(N-1):-0.5];
X = abs(X).*2;
cneg = 1-[0 0 1];
cpos = 1-[1 0 0];
R = [1-X(:,1:N/2)*cneg(1) 1-X(:,(N/2 + 1):N)*cpos(1)];
G = [1-X(:,1:N/2)*cneg(2) 1-X(:,(N/2 + 1):N)*cpos(2)];
B = [1-X(:,1:N/2)*cneg(3) 1-X(:,(N/2 + 1):N)*cpos(3)];
cmap = [R' G' B'];

cmap2 = cmap;

%% make iaMat from the edges
edges = subGraphStruct.subGraphDetails(idx).edges;
iaMat = zeros(max(max(edges)));
for i = 1:size(edges)
    iaMat(edges(i,1),edges(i,2)) = 1;
end
%redIaMat

uCmpdIdx = subGraphStruct.subGraphDetails(idx).cmpdIdx;
iaMat = iaMat(uCmpdIdx,:);
iaMat = iaMat(:,uCmpdIdx);
id2collapse = [];
cmpdName = {''};
for i = 1:length(uCmpdIdx)
    if mns.observations.isMet(uCmpdIdx(i))
        cmpdName(i) = metName(uCmpdIdx(i)-mns.observations.noOfRp);
        nodeId(i) = metName(uCmpdIdx(i)-mns.observations.noOfRp);
        shape{i} = 'ellipse';
    else
        temp = mns.observations.iaId2dataId{uCmpdIdx(i)};
        nodeId(i) = rpName(uCmpdIdx(i));
        tempName = unique(geneName(temp));
        tempStr = tempName{1};
        for j = 2:length(tempName)
            tempStr = [tempStr '; ' tempName{j}];
        end
        
        if collapseGeneIds
            idExist = find(strcmp(cmpdName, tempStr) == 1);
            if ~isempty(idExist)
                iaMat(idExist,:) = iaMat(idExist,:)+iaMat(i,:);
                iaMat(:,idExist) = iaMat(:,idExist)+iaMat(:,i);
                id2collapse = [id2collapse,i];
                cmpdName(i) = {''};
            else
                cmpdName(i) = {tempStr};
            end
        else
            cmpdName(i) = {tempStr};
        end
        
        shape{i} = 'box';
    end
end
%% if collapse gene ids
if collapseGeneIds
    iaMat(iaMat > 1) = 1;
    iaMat(id2collapse,:) = [];
    iaMat(:,id2collapse) = [];
    cmpdName(id2collapse) = [];
    nodeId(id2collapse) = [];
    shape(id2collapse) = [];
    uCmpdIdx(id2collapse) = [];
end

%% make biograph object
BGObj2 = biograph(iaMat, nodeId,'LayoutType','equilibrium');
%% make Size and color vlaues
pval = zeros(length(uCmpdIdx),1);
fc = zeros(length(uCmpdIdx),1);
detected = mns.observations.detected(uCmpdIdx);
for i = 1:length(uCmpdIdx)
    if detected(i);
        if mns.observations.isMet(uCmpdIdx(i))
            idxTemp = mns.observations.iaId2dataId{uCmpdIdx(i)};
            pval(i) = pvalData(idxTemp);
            fc(i) = fcData(idxTemp);
        else 
            pval(i) = mns.observations.data(uCmpdIdx(i));
            tempFc = fcGeneData(mns.observations.iaId2dataId{uCmpdIdx(i)});
            fc(i) = median(tempFc);
            detected(i) = 1;
        end
    end
end
Min_all_pval = 0;
Max_all_pval = 10;
cmpdSizeMax = 10;
cmpdSizeMin = 2;
cmpdSizeValTruncate = 10;
sizeVals = -log10(pval);
sizeVals(isnan(sizeVals)) = min(sizeVals) - abs(min(sizeVals));
sizeVals(isinf(sizeVals)) = NaN;
sizeVals(isnan(sizeVals)) = max(sizeVals) + abs(max(sizeVals));
sizeVals(abs(sizeVals)>abs(cmpdSizeValTruncate)) = cmpdSizeValTruncate;%*sign(sizeVals(abs(sizeVals)>abs(cmpdSizeValTruncate)));
cmpdSize = round((sizeVals-Min_all_pval)/(Max_all_pval-Min_all_pval)*(cmpdSizeMax-cmpdSizeMin)+cmpdSizeMin);
cmpdSize(~isfinite(cmpdSize)) = cmpdSizeMin;
%%
% cmap = redgreencmap(255);
colorVals = fc;
colorVals(isnan(colorVals)) = min(colorVals) - abs(min(colorVals));
colorVals(isinf(colorVals)) = NaN;
colorVals(isnan(colorVals)) = max(colorVals) + abs(max(colorVals));
maxVal = cmax;
cmpdColorAxisMin = -maxVal;
cmpdColorAxisMax = maxVal;
colors = round((colorVals-cmpdColorAxisMin)/(cmpdColorAxisMax-cmpdColorAxisMin)*253)+1;
colors(colors<1)= 1;
colors(colors>254)= 254;
cmpdColor = cmap(colors,:);
%%

for i = 1:length(uCmpdIdx)
    if detected(i) == 1
        BGObj2.Nodes(i).Color = cmpdColor(i,:);
    else
        BGObj2.Nodes(i).Color = [0.5 0.5 0.5];
    end
    BGObj2.Nodes(i).LineColor = [0.3 0.3 0.3];
    BGObj2.Nodes(i).Size = [cmpdSize(i) cmpdSize(i)];
    BGObj2.Nodes(i).Label = cmpdName{i};
    BGObj2.Nodes(i).Shape = shape{i};
end

view(BGObj2)

end