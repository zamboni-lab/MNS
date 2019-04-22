function mns_initLabels2labelFile(initLabels, filename)

% open outfile for observable edge data
fid = fopen(filename, 'w');
fprintf(fid,'%d\r\n', length(initLabels), length(groupMean));
c = 0;
for i = 1:length(metIdIA)
    met.modelVarId(i) = i-1;
    idx = find(aLink2Ia(:,1) == i);
    if ~isempty(idx)
        c = c+1;
        met.obsFacId(i) = c;
        met.detected(i) = 1;
        idTemp = dataId(aLink2Ia(idx,2));
        met.iaId2dataId(i) = idTemp(1); %shoul be as well checked with score
        dataTemp = data(dataId(aLink2Ia(idx,2)));
        met.data(i) = dataTemp(1); %Should be compared with score
        met.obsCliquePotential(i,:) = mns_calcProbGaussian(groupMean, groupStd, met.data(i),obsFuncType, l2, nCliques.nCliquePerVar);
    end
    if p2file
        for j = -1:length(groupMean)
            if j == -1
                fprintf(fid,'%d', met.modelVarId(i));
            elseif j == 0
                fprintf(fid,'\t%d', met.detected(i));
            else
                fprintf(fid,'\t%f', met.obsCliquePotential(i,j));
            end
            
        end
        
        if i ~= length(metIdIA)
            fprintf(fid,'\r\n');
        end
    end
end

if p2file
    fclose(fid);
end
end