function out = mns_subGraphMetGeneResult2table(subGraphResults, model, uIonsCo, rpoFile)

if nargin < 3
    uIonsCo = 0;
end
if nargin < 4
    doRpoPtwAnalysis = false;
    ptwHead = [];
    ptwDetails = [];
else
    doRpoPtwAnalysis = true;
    ptwHead = {'Pathway Contribution'};
    [ pws, rps ] = parseOntology( rpoFile );
    pwDesc = regexprep(regexprep([{pws(:).description};]', ' - Homo .*',''), '^\s','');
end

noUIonsCol = 16;

subGraphHead = subGraphResults.mergedSubGraphHeader;
subGraphList = subGraphResults.mergedSubGraphList;

%make List of metabolite ids and names for each
%subGraph
metName = cell(size(subGraphList,1),1);
geneName = cell(size(subGraphList,1),1);
metId = cell(size(subGraphList,1),1);
% col2parse = [5,6,7,9,12];
% col2parse = [8, 9, 10, 11, 12, 14, 17];
ptwCo = 0.1;
%% 
noOfGeneCmpds = max(find(model.isRp == 1));
rowIdx = [];
for j = 1:length(metName)
    if cell2mat(subGraphList(j,16)) < uIonsCo
        continue;
    end
    rowIdx = [rowIdx; j];
    subMetIdx = subGraphResults.mergedSubGraphStruct.subGraphDetails(j).metIdx;
    subRpIdx = subGraphResults.mergedSubGraphStruct.subGraphDetails(j).rpIdx;
    rpId = regexprep(model.rpId(subRpIdx), '\_\d+','');
    for k = 1:length(subMetIdx)
        if k == 1
            metName{j} = model.metaboliteName{subMetIdx(k)-noOfGeneCmpds};
            metId{j} = model.metaboliteId{subMetIdx(k)-noOfGeneCmpds};
        else
            metName{j} = [metName{j} ', ' model.metaboliteName{subMetIdx(k)-noOfGeneCmpds}];
            metId{j} = [metId{j} ', ' model.metaboliteId{subMetIdx(k)-noOfGeneCmpds}];
        end
    end
    if doRpoPtwAnalysis
        ptwTemp = zeros(length(pws),1);
        for i = 1:length(pws)
            ptwTemp(i) = length(intersect(rps{i},rpId));
        end
        ptwTemp = ptwTemp/sum(ptwTemp);
        idxPtw = find(ptwTemp >= ptwCo);
        ptwDetails{j,1} = '';
        if ~isempty(idxPtw)
            [~,idxS] = sort(ptwTemp(idxPtw), 'descend');
            for i = 1:length(idxPtw)
                ptwDetails{j,1} = [ptwDetails{j,1} sprintf('%s (%d%%); ', pwDesc{idxPtw(idxS(i))}, round(ptwTemp(idxPtw(idxS(i)))*100))];
            end
        end
        
    end
    geneSymbTemp = [];
    for k = 1:length(subRpIdx)
        geneSymbTemp = [geneSymbTemp model.GeneSymbol(model.gene2rp(:,subRpIdx(k)) == 1)];
    end
    geneSymbTemp = unique(geneSymbTemp);
    for k = 1:length(geneSymbTemp)
        if k == 1
            geneName{j} = geneSymbTemp{k};
%             metId{j} = model.metaboliteId{subMetIdx(k)-noOfGeneCmpds};
        else
            geneName{j} = [geneName{j} ', ' geneSymbTemp{k}];
%             metId{j} = [metId{j} ', ' model.metaboliteId{subMetIdx(k)-noOfGeneCmpds}];
        end
    end
%     for k = 1:length(col2parse)
%         parseData = subGraphList{j,col2parse(k)};
%         for l = 1:size(parseData,1);
%             if l == 1
%                 subGraphList{j,col2parse(k)} = num2str(parseData(l));
%             else
%                 subGraphList{j,col2parse(k)} = [subGraphList{j,col2parse(k)} ', ' num2str(parseData(l))];
%             end
%         end
%     end
end
%%
col2export1 = [1 2 3 4 5 6 7];
col2export2 = [13 15 16];
outHeader = [subGraphHead(col2export1) {'MetId'} {'MetName'} {'GeneName'} subGraphHead(col2export2) ptwHead];
if doRpoPtwAnalysis
    ptwDetails = ptwDetails(rowIdx);
end
outData = [subGraphList(rowIdx,col2export1) metId(rowIdx) metName(rowIdx) geneName(rowIdx)...
        subGraphList(rowIdx,col2export2) ptwDetails];

out = [outHeader; outData];

end


function [ pws, cmpds ] = parseOntology( filename )
%PARSEMO Summary of this function goes here
%   Detailed explanation goes here
fh = fopen(filename);
pwCount = 0;
cmpds = {};
curr_cmpds = [];
while 1
    tline = fgetl(fh);
    if isnumeric(tline) || isempty(tline)
        if ~isempty(curr_cmpds)
            cmpds(pwCount) = {curr_cmpds};
        else
            cmpds(pwCount) = {' '};
        end
        break;
    end
    switch tline(1)
        case '>' % new pathway
            if ~isempty(curr_cmpds)
                cmpds(pwCount) = {curr_cmpds};
            end
            pwCount = pwCount + 1;
            curr_cmpds = [];
            pos = strfind(tline, '#');
            pws(pwCount).ID = tline(2:pos-1);
            pws(pwCount).description = tline(pos+1:end);
        case 'c' % compound
            curr_cmpds = [curr_cmpds {tline(5:end)}];
        case 'r' % compound
            switch tline(1:2)
                case 'rp'
                    curr_cmpds = [curr_cmpds {tline(4:end)}];
                case 'rc'
                    curr_cmpds = [curr_cmpds {tline(6:end)}];
            end
    end
end
fclose(fh);
end