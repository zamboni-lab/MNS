function nFracListOut = mns_filterSeqResultsNFracList(nFracList, minNoFractures, minRatioFracturesToGroups)

if nargin < 2
    minNoFractures = 2;
end
if nargin < 3
    minRatioFracturesToGroups = 2;
end

noFrac = cell2mat(nFracList(2:end,8));
fracGroups = cell2mat(nFracList(2:end,9));
ratio = noFrac./fracGroups;

idx = find(ratio >= minRatioFracturesToGroups & noFrac >= minNoFractures);

nFracListOut = nFracList([1; idx+1],:);
nFracListOut(2:end,1) = num2cell(tiedrank(-cell2mat(nFracListOut(2:end,10))));
end