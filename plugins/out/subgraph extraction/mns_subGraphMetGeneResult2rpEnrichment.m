function out = mns_subGraphMetGeneResult2rpEnrichment(subGraphResults, model, uIonsCo, rpoFile)

%% parse rpo file
[ pws, rps ] = parseOntology( rpoFile );
pwNames = regexprep(regexprep([{pws(:).description};]', ' - Homo .*',''), '^\s','');

%% total numer or rp
rp = regexprep(model.rpId, '\_\d+','');
uRp = unique(rp);
nRpTotal = length(uRp);
%% get subGraphList
subGraphList = subGraphResults.mergedSubGraphList;

idxSubGraphs = find(cell2mat(subGraphList(:,16)) >= uIonsCo);

%% calculate enrichment
enrichScore = ones(length(pwNames), length(idxSubGraphs));

for i = 1:length(pwNames)
    nRpPtw = length(rps{i});
    for j = 1:length(idxSubGraphs)
        subRpIdx = subGraphResults.mergedSubGraphStruct.subGraphDetails(idxSubGraphs(j)).rpIdx;
        uRpId = unique(regexprep(model.rpId(subRpIdx), '\_\d+',''));
        nRpSubGraph = length(uRpId);
        nRpIntersect = length(intersect(uRpId, rps{i}));
        
        tmp = hygepdf(1:nRpPtw, nRpTotal, nRpPtw, nRpSubGraph);
        score = hygepdf(nRpIntersect, nRpTotal, nRpPtw, nRpSubGraph);
        apex = find(tmp==max(tmp));
        if nRpIntersect>=apex
            if score < enrichScore(i,j)
                enrichScore(i,j) = score;
            end
        end
    end
end
%% correct for multiple testing
for c = 1:size(enrichScore,2)
    try
        [enrichmentScoresStoreyFDR(:,c),enrichmentScoresStorey(:,c)] = mafdr(enrichScore(:,c));
    catch
        enrichmentScoresStoreyFDR(:,c) = NaN;
        enrichmentScoresStorey(:,c) = NaN;
    end
    
    [~,~,enrichmentScoresAdjPvalBH(:,c)] = fdr(enrichScore(:,c));
end
%% generate outStructure
out.enrichmentScores = enrichScore;
out.enrichmentScoresQvalues = enrichmentScoresStorey;
out.pathways = pwNames;
out.pairs = num2cell(idxSubGraphs);
out.overviewQ = [{''} out.pairs'; pwNames num2cell(enrichmentScoresStorey)];
out.overview = [{''} out.pairs'; pwNames num2cell(enrichScore)];
% out.overviewSigHitsSinglePtw = overviewSigHitsSinglePtw;
% out.overviewSigHitsTotal = overviewSigHitsTotal;
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