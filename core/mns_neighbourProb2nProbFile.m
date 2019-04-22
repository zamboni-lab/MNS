%% Copyright and License
% Metabolic Network Segmentation Toolbox – A probabilistic graphical 
% modelling tool to identify sites and sequential order of metabolic regulations
% Copyright (C) 2016, Andreas Kühne & Nicola Zamboni
% 
% This file is part of the Metabolic network segmentation toolbox
% 
% Metabolic network segmentation toolbox is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation, either
% version 3 of the License or any later version.
% 
% Metabolic network segmentation toolbox is distributed in the hope 
% that it will be useful, but WITHOUT ANY WARRANTY; without even the
% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy (gpl.txt) of the GNU General Public License
% along with the metabolic network segmentation toolbox. If not,
% see http://www.gnu.org/licenses/.
%
function out = mns_neighbourProb2nProbFile(filename, l1, minClSize, maxClSize, nFuncType, p2file)
%nFuncType
% 0 -> weighted;
% 1 -> weighted + highest score independent of clique size;
% 2 -> not weighted;
% 3 -> not weighted + highest score independent of clique size;
    
if nargin < 5
    weighted = 0;
end
prob = nan(maxClSize-minClSize+1,maxClSize);
uLabels = nan(maxClSize-minClSize+1,maxClSize);
nCl = length(minClSize:maxClSize);
if p2file
    fid = fopen(filename, 'w');
end
% Format
% 1st Line: No CliqueTemplates, minCliqueSize, maxCliqueSize
% followed: cliqueSize, total perm, all different, 1 different ..., all
% same
% fidPerm = fopen(filenamePerm, 'w');

if p2file
    fprintf(fid, '%d\t%d\t%d\r\n',nCl,minClSize,maxClSize);
end
%fprintf(fidPerm, '%d\t%d\t%d\r\n',nCl,minClSize,maxClSize);


for i = minClSize:maxClSize
    switch nFuncType
        case 0
            [prob(i-minClSize+1, 1:i), uLabels(i-minClSize+1, 1:i)] = calcExpProb(i,l1, false);
        case 1
            [prob(i-minClSize+1, 1:i), uLabels(i-minClSize+1, 1:i)] = calcExpProb(i,l1, true);
        case 2
            [prob(i-minClSize+1, 1:i), uLabels(i-minClSize+1, 1:i)] = calcExpProbNotWeighted(i,l1, false);
        case 3
            [prob(i-minClSize+1, 1:i), uLabels(i-minClSize+1, 1:i)] = calcExpProbNotWeighted(i,l1, true);
    end
    if p2file
        fprintf(fid,'%d',i);
    end
    out.cliqueTemplateSize(i-minClSize+1) = i;
%     prob = prob*8;
    for j = 1:i
        if p2file
            fprintf(fid,'\t%f',prob(i-minClSize+1, j));
        end
        if j == 1
            n(j) = 1;
        elseif j == i
            n(j) = 1;
        else
            n(j) = calcNoPerm(j, i);
        end
    end
%     fprintf(fidPerm,'%d\t%d',i,sum(n));
    for j = 1:i
%         fprintf(fidPerm,'\t%d',n(j));
    end
    
    clear n
    if p2file
        if i ~= maxClSize
            fprintf(fid,'\r\n');
            %         fprintf(fidPerm,'\r\n');
        end
    end
end
if p2file
    fclose(fid);
end
out.uniqueLabels = nanmax(uLabels,[],1);
out.noOfcliqueTemplates = nCl;
out.probMat = prob;
out.uniqueLabelsMat = uLabels;


end

function n = calcNoPerm(picked, setSize)
n = factorial(setSize)/( factorial(setSize-picked)*factorial(picked) );
end

function [prob2, uL] = calcExpProb(clSize, l1, cliqueIndep)
prob2 = zeros(clSize,1);
uL = zeros(clSize,1);
uCl = 1:clSize;

if ~cliqueIndep
    cld = 0;
else
    cld = 1;
end

for i = 1:clSize
    uL(i) = i;
    prob2(i) = exp(-l1*((uCl(i)-cld)/clSize));
end

end

function [prob2, uL] = calcExpProbNotWeighted(clSize, l1, cliqueDep)
prob2 = zeros(clSize,1);
uL = zeros(clSize,1);
uCl = 1:clSize;
if cliqueDep
    cld = 0;
else
    cld = 1;
end
for i = 1:clSize
    uL(i) = i;
    prob2(i) = exp(-l1*(uCl(i)-cld));
end

end