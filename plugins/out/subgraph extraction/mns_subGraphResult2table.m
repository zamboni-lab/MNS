function out = mns_subGraphResult2table(subGraphResults, model)

subGraphHead = subGraphResults.mergedSubGraphHeader;
subGraphList = subGraphResults.mergedSubGraphList;

%make List of metabolite ids and names for each
%subGraph
metName = cell(size(subGraphList,1),1);
metId = cell(size(subGraphList,1),1);
col2parse = [5,6,7,9,12];
for j = 1:length(metName)
    subMetIdx = subGraphResults.mergedSubGraphStruct.subGraphDetails(j).metIdx;
    for k = 1:length(subMetIdx)
        if k == 1
            metName{j} = model.metaboliteName{subMetIdx(k)};
            metId{j} = model.metaboliteId{subMetIdx(k)};
        else
            metName{j} = [metName{j} ', ' model.metaboliteName{subMetIdx(k)}];
            metId{j} = [metId{j} ', ' model.metaboliteId{subMetIdx(k)}];
        end
    end
    for k = 1:length(col2parse)
        parseData = subGraphList{j,col2parse(k)};
        for l = 1:size(parseData,1);
            if l == 1
                subGraphList{j,col2parse(k)} = num2str(parseData(l));
            else
                subGraphList{j,col2parse(k)} = [subGraphList{j,col2parse(k)} ', ' num2str(parseData(l))];
            end
        end
    end
end
outHeader = [subGraphHead(1:8) {'MetModelId'} {'MetId'} {'MetName'} subGraphHead(10:end)];
outData = [subGraphList(:,1:9) metId metName subGraphList(:,10:end)];

out = [outHeader; outData];

end