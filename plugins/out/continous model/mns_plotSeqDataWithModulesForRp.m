function mns_plotSeqDataWithModulesForRp(mnsResults, mnsDataStruct, model, idx, rpIdx, plotMode)
if nargin < 6 || isempty(plotMode)
    plotMode = 'area';
end

tempResults = mnsResults.mns{1};
nTp = size(mnsDataStruct.data,2);

for i = 1:length(rpIdx)
    metIdx = find(model.mat(rpIdx(i),:) ~= 0);
    mns_plotTcWithModules(mnsResults, mnsDataStruct, model, idx, plotMode, metIdx)
end

end