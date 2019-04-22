function modelMNS = mns_tbModel2MnsModel(model)

oldNames = {'Compounds' 'CompoundsName' 'RP' 'CompoundstoCompoundsviaRP' 'CompoundstoRP'};
newNames = {'metaboliteId', 'metaboliteName' 'rpId' 'iaMat' 'mat'};

modelMNS = model;
for i = 1:length(oldNames)
    modelMNS.(newNames{i}) = modelMNS.(oldNames{i});
    rmfield(modelMNS, oldNames{i});
end

modelMNS.mat = modelMNS.mat';
[modelMNS.gene2geneDist,modelMNS.gene2rp] = mns_gene2geneDistMat(modelMNS);
end
