function out = mns_subGraphMetGeneResult2frequencyPlots(subGraphResults, model, uIonsCo, rpoFile)

%% parse rpo file
[ pws, rps ] = parseOntology( rpoFile );
pwNames = regexprep(regexprep([{pws(:).description};]', ' - Homo .*',''), '^\s','');
%% get subGraphList
subGraphList = subGraphResults.mergedSubGraphList;
idxSubGraphs = find(cell2mat(subGraphList(:,16)) >= uIonsCo);

%% go through all subgraphs and count met, genes, and edges

noOfMet = cell2mat(subGraphList(idxSubGraphs,6));
noOfGenes = cell2mat(subGraphList(idxSubGraphs,13));
noOfEdges = cell2mat(subGraphList(idxSubGraphs,2));
ptwFreq = zeros(length(pws), length(idxSubGraphs));
for j = 1:length(idxSubGraphs)
%     rowIdx = [rowIdx; j];
%     subMetIdx = subGraphResults.mergedSubGraphStruct.subGraphDetails(idxSubGraphs(j)).metIdx;
%     noOfMet(j) = length(unique(subMetIdx));
%     noOfEdges(j) = length(unique(subMetIdx));
    subRpIdx = subGraphResults.mergedSubGraphStruct.subGraphDetails(idxSubGraphs(j)).rpIdx;
    rpId = regexprep(model.rpId(subRpIdx), '\_\d+','');
%     ptwTemp = zeros(length(pws),1);
    for i = 1:length(pws)
        ptwFreq(i,j) = length(intersect(rps{i},rpId));
    end
    ptwFreq(:,j) = ptwFreq(:,j)./sum(ptwFreq(:,j));
        

end

%% define out
out.pathways = pwNames;
out.pairs = num2cell(idxSubGraphs);
out.noOfMets = noOfMet;
out.noOfGenes = noOfGenes;
out.noOfEdges = noOfEdges;
out.ptwFreq = ptwFreq;

%% plot the data
cAxisMax = 20;
figure,
subplot(3,1,1);
imagesc(out.noOfMets')
ylabel('#Metabolites')
set(gca,'XTick', 1:length(idxSubGraphs), 'XTickLabel', idxSubGraphs);
set(gca,'YTick', []);
caxis([0 cAxisMax]);
subplot(3,1,2);
imagesc(out.noOfGenes')
ylabel('#Genes')
set(gca,'XTick', 1:length(idxSubGraphs), 'XTickLabel', idxSubGraphs);
set(gca,'YTick', []);
caxis([0 cAxisMax]);
subplot(3,1,3);
imagesc(out.noOfEdges')
ylabel('#edges')
set(gca,'XTick', 1:length(idxSubGraphs), 'XTickLabel', idxSubGraphs);
set(gca,'YTick', []);
caxis([0 cAxisMax]);
%% make colormap
offset = 0;
upOffset = 0;

cmax = cAxisMax;
cSize = cmax;
cpos = [0 102 204]/255;
cpos = 1-cpos;
cmap = ones(cSize,3);
cmap(offset+1:cSize,1) = 1-[0:1/(cSize-offset-1)*cpos(1):1*cpos(1)];
cmap(offset+1:cSize,2) = 1-[0:1/(cSize-offset-1)*cpos(2):1*cpos(2)];
cmap(offset+1:cSize,3) = 1-[0:1/(cSize-offset-1)*cpos(3):1*cpos(3)];
if upOffset > 0
   cmap(cSize+1:cSize+offset,1) = 1; 
end
colormap(cmap)

h = colorbar;
% set(get(h,'YLabel'),'String','-log_{10}(q-value)')

clear h offset upOffset cmap cSize  cmax
%% plot ptwFrequency
figure
imagesc(out.ptwFreq*100)
set(gca,'XTick', 1:length(idxSubGraphs), 'XTickLabel', idxSubGraphs);
set(gca, 'YTick', 1:length(pwNames), 'YTickLabel',  pwNames);

%% make colormap for frequency
offset = 0;
upOffset = 0;

cmax = 100;
cSize = cmax;
cpos = [204 0 51]/255;
cpos = 1-cpos;
cmap = ones(cSize,3);
cmap(offset+1:cSize,1) = 1-[0:1/(cSize-offset-1)*cpos(1):1*cpos(1)];
cmap(offset+1:cSize,2) = 1-[0:1/(cSize-offset-1)*cpos(2):1*cpos(2)];
cmap(offset+1:cSize,3) = 1-[0:1/(cSize-offset-1)*cpos(3):1*cpos(3)];
if upOffset > 0
   cmap(cSize+1:cSize+offset,1) = 1; 
end
colormap(cmap)

h = colorbar;
set(get(h,'YLabel'),'String','Reaction pairs in pathway')
caxis([0 cmax]);
% clear h offset upOffset cmap cSize cmax

%% find ptw not covered
freqCo = 0.1;
temp = zeros(size(out.ptwFreq));
temp(out.ptwFreq >= freqCo) = 1;
idxPtw = find(sum(temp,2) ~= 0);
figure
imagesc(out.ptwFreq(idxPtw,:)*100)
set(gca,'XTick', 1:length(idxSubGraphs), 'XTickLabel', idxSubGraphs);
set(gca, 'YTick', 1:length(pwNames(idxPtw)), 'YTickLabel',  pwNames(idxPtw));
h = colorbar;
set(get(h,'YLabel'),'String','Reaction pairs in pathway')
caxis([0 cmax]);
colormap(cmap)
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