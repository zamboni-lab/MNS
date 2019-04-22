function [modelOut, modIdx, feIdx ] = mns_reducedModelAndDataByPtw(ptwStruct, ptwId, model, fe)
% Input variables
% ptwStruct: struct with format ptwStruct(ptwId).ptw.cmpdId
% ptwId: int array with ids of selected pathways
% model: metabolic model including iaMat and metabolite Id
% fe: fiaExp


% get unique metabolite ids of pathways
pIdx = []; %cmpd Ids in a given Ptw
for i = 1:length(ptwId)
    pIdx = [pIdx;ptwStruct(ptwId(i)).ptw.cmpdID];
end

pIdx = unique(pIdx);

% make annotation map
aId = [];
aScore = [];
aName = [];
aFeIdx = [];

for i = 1:length(fe.merged.annotation)
    aId = [aId; fe.merged.annotation(i).id'];
    aScore = [aScore; fe.merged.annotation(i).score'];
    aName = [aName; fe.merged.annotation(i).name'];
    aFeIdx = [aFeIdx; repmat(i,fe.merged.annotation(i).nHits,1)];
end


%map pIdx to model and fiaExp indices
modIdx = zeros(length(pIdx),1);
feIdx = [];
for i = 1:length(pIdx)
    modIdx(i) = find(strcmp(model.metaboliteId, pIdx{i}) == 1);
    tempFeIdx = find(strcmp(aId, pIdx{i}) == 1);
    if ~isempty(tempFeIdx)
        feIdx = [feIdx; aFeIdx(tempFeIdx(aScore(tempFeIdx) == max(aScore(tempFeIdx))))];
    end
end


% reduce Model
modelOut.metaboliteId = pIdx;
modelOut.iaMat = model.iaMat(modIdx, modIdx);

if isfield(model, 'metaboliteName')
    modelOut.metaboliteName = model.metaboliteName(modIdx);
end
if isfield(model, 'metaboliteFormula')
    modelOut.metaboliteFormula = model.metaboliteFormula(modIdx);
end
if isfield(model, 'metaboliteExactMass')
    modelOut.metaboliteExactMass = model.metaboliteExactMass(modIdx);
end
if isfield(model, 'metaboliteMolWeight')
    modelOut.metaboliteMolWeight = model.metaboliteMolWeight(modIdx);
end





end