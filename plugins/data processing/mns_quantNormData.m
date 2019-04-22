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
function dataOut = mns_quantNormData(data, qStep)

idxNeg = find(data < 0);
idxPos = find(data > 0);

quant = 0:qStep:1;
quantPos = zeros(length(quant),1);
quantNeg = zeros(length(quant),1);
quantMean = quantPos;

for i = 1:length(quant)
    quantPos(i) = quantile(data(idxPos), quant(i));
    quantNeg(i) = quantile(abs(data(idxNeg)), quant(i));
    quantMean(i) = mean([quantPos(i), abs(quantNeg(i))]);
end
%% do linear transform within quantiles
dataOut = data;
for i = 1:length(idxPos)
    idxLow = find(quantPos <= data(idxPos(i)));
    idxHigh = find(quantPos >= data(idxPos(i)));
    if idxLow(end) == idxHigh(1)
        dataOut(idxPos(i)) = quantMean(idxLow(end));
    else
        m = (quantMean(idxHigh(1))-quantMean(idxLow(end)))/(quantPos(idxHigh(1))-quantPos(idxLow(end)));
        c = quantMean(idxLow(end))-m* quantPos(idxLow(end));
        dataOut(idxPos(i)) = data(idxPos(i))*m+c;
    end
end

for i = 1:length(idxNeg)
    idxLow = find(quantNeg <= abs(data(idxNeg(i))));
    idxHigh = find(quantNeg >= abs(data(idxNeg(i))));
    if idxLow(end) == idxHigh(1)
        dataOut(idxNeg(i)) = -quantMean(idxLow(end));
    else
        m = (quantMean(idxHigh(1))-quantMean(idxLow(end)))/(quantPos(idxHigh(1))-quantPos(idxLow(end)));
        c = quantMean(idxLow(end))-m* quantPos(idxLow(end));
        dataOut(idxNeg(i)) = -(data(idxNeg(i))*m+c);
    end
end

end