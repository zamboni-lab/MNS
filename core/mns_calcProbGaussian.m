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
function [prob] = mns_calcProbGaussian(meanVec, stdVec, data, obsFuncType,l2,nCliquePerVar,neighborMetabolites, l1)

% obsFuncType
% 1: %normalized gaussian
% 2: %gaussian
% 3: %normalized and linkage dependendt gaussian
% 4: %linkage dependendt gaussian


prob = zeros(length(meanVec),1);
meanMax = max(meanVec);
meanMin = min(meanVec);
for i = 1:length(prob)
    switch obsFuncType
        case 1 %normalized gaussian
            prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
        case 2 %gaussian
            prob(i) = exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
        case 3 %normalized and linkage dependendt gaussian
            prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( -l2*nCliquePerVar* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
        case 4 %linkage dependendt gaussian
            prob(i) = exp( -l2 * nCliquePerVar * (data-meanVec(i)).^2 / (2*stdVec(i)^2));
        case 5 %gaussian with borders = 1
            prob(i) = exp( -l2* (data-meanVec(i)).^2 / (2*stdVec(i)^2));
            if meanVec(i) == meanMin
                if data <= meanMin
                    prob(i) = 1;
                end
            elseif meanVec(i) == meanMax
                if data >= meanMax
                    prob(i) = 1;
                end
            end
        case 6 %neighbor metabolite dependent and normalized gaussian
            prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * exp( -(1+(neighborMetabolites)*l2)*(data-meanVec(i)).^2 / (2*stdVec(i)^2));
        case 7 %neighbor metabolite dependent  gaussian
%             prob(i) = exp( -(1+l2*neighborMetabolites*l1) * (data-meanVec(i)).^2 / (2*stdVec(i)^2));
            prob(i) = exp( -(1+l2*2^(neighborMetabolites-1)) * (data-meanVec(i)).^2 / (2*stdVec(i)^2));
		case 8 %neighbor metabolite dependent II and normalized gaussian
            prob(i) = 1/(sqrt(2*pi)* stdVec(i) ) * l2*neighborMetabolites * exp( -l2*(data-meanVec(i)).^2 / (2*stdVec(i)^2));
		case 9 %neighbor metabolite dependent II and normalized gaussian
            prob(i) = l2*(1+neighborMetabolites-1) * exp( -l2*(data-meanVec(i)).^2 / (2*stdVec(i)^2));
    end
    if isnan(prob(i))
        prob(i) = 0;
    end
end

end