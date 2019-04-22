function mns_subGraphResult2bioGraph(subGraphStruct, pvalData,fcData,mns,metName,idx, cmax)
if nargin < 6
    idx = 1;
end
if nargin < 7
    cmax = 2.5;
end

%% make iaMat from the edges
edges = subGraphStruct.subGraphDetails(idx).edges;
iaMat = zeros(max(max(edges)));
for i = 1:size(edges)
    iaMat(edges(i,1),edges(i,2)) = 1;
end
%redIaMat
uMetIdx = subGraphStruct.subGraphDetails(idx).metIdx;
iaMat = iaMat(uMetIdx,:);
iaMat = iaMat(:,uMetIdx);
metName = metName(uMetIdx);
%% make biograph object
BGObj2 = biograph(iaMat, metName,'LayoutType','equilibrium');
%% make Size and color vlaues
dataIdx = mns.observations.iaId2dataId(uMetIdx);
pval = zeros(length(uMetIdx),1);
fc = zeros(length(uMetIdx),1);
detected = zeros(length(uMetIdx),1);
for i = 1:length(dataIdx)
    if dataIdx(i) ~= 0
        pval(i) = pvalData(dataIdx(i));
        fc(i) = fcData(dataIdx(i));
        detected(i) = 1;
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
%% make dualColorColormap
N = 255;
X = [0.5: -1/(N-1):-0.5];
X = abs(X).*2;
cneg = 1-[0 1 0];
cpos = 1-[1 0 0];
R = [1-X(:,1:N/2)*cneg(1) 1-X(:,(N/2 + 1):N)*cpos(1)];
G = [1-X(:,1:N/2)*cneg(2) 1-X(:,(N/2 + 1):N)*cpos(2)];
B = [1-X(:,1:N/2)*cneg(3) 1-X(:,(N/2 + 1):N)*cpos(3)];
cmap = [R' G' B'];

cneg = [0 102 204]/255;
cpos = [204 0 51]/255;

n = 255;

cmap = exp_colormap('custom white', n, cneg, cpos);
% figure

%%
% cmax = 2.5;
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
BGObj2.ShowArrows = 'off';
for i = 1:length(uMetIdx)
    if detected(i) == 1
        BGObj2.Nodes(i).Color = cmpdColor(i,:);
        BGObj2.Nodes(i).LineColor = [0.2 0.2 0.2];
    else
        BGObj2.Nodes(i).Color = [0.5 0.5 0.5];
        BGObj2.Nodes(i).LineColor = [0.2 0.2 0.2];
    end
    BGObj2.Nodes(i).Size = [cmpdSize(i) cmpdSize(i)];
end

view(BGObj2)

end