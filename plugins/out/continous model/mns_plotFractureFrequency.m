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
function mns_plotFractureFrequency(mnsResults, model,wsMax,wnMax,dWs, dWn)

wsVec = 0:dWs:wsMax;
wnVec = 0:dWn:wnMax;

%%
sFractures = length(wnVec);
nFractures = length(wnVec);
scoreVec = length(wnVec);
probVec = length(wnVec);
for i = 1:length(wnVec)
    [prob, score, idx,idxScore] = ...
        mns_calcProbabilityScanData(mnsResults, model,wsVec(i),wnVec(i),false,false);
    %         mns_plotTemporalScreeningResults(mnsResults,prob, score,true)
    idxRow = find(mnsResults.nL1Range == mnsResults.nL1(idx));
    idxCol = find(mnsResults.tL3Range == mnsResults.tL3(idx));
    nFractures(i) = mnsResults.amountBreakingReaction(idxRow, idxCol);
    sFractures(i) = mnsResults.amountTemporalBreaks(idxRow, idxCol);
    scoreVec(i) = score(idxScore(1),idxScore(2));
    probVec(i) = prob(idxScore(1),idxScore(2));
end
%%
figure
subplot(1,2,1)
hold all
plot(wsVec, nFractures, '-')
plot(wsVec, sFractures, '-')
legend('neighborhood fractures','sequence fractures')
ylabel('#fractures');
xlabel('ws, wn');
% figure
% hold all
% plot(wtVec, scoreVec, '-')
% ylabel('Score');
% xlabel('ws, wn');
subplot(1,2,2)
hold all
plot(wsVec, probVec, '-')
ylabel('Sum of Observation potential');
xlabel('ws, wn');

end