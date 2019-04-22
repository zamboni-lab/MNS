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
function initMNS = mns_generateInitMNS(varargin)

initMNS.labeling = 'kmeans';
initMNS.noOfclusters = 3;
% initMNS.noOfclusters = 5; previously 5
initMNS.weights = 'exact';
initMNS.oL2max = 1;
initMNS.oL2min = 1;
initMNS.nL1max = 1;
initMNS.nL1min = 1;
initMNS.nL1steps = 40;
initMNS.tL3min = 1; % if exact -> oL2min = oL2; if gaussian -> oL2min = mean
initMNS.tL3max = 1; % if exact -> oL2max = oL2; if gaussian -> oL2max =
initMNS.tL3steps = 40;
% initMNS.meanType = 'initLabels'; previous
initMNS.meanType = 'initLabels';
%options:
% 'initLabels' - mean determined by groups defined by the initial labels e.g. through k-means
% 'linear - std' - means linearly distributed between min(data)+std and max(data)-std
% 'linear - quantile' - means linearly distributed between quantile(data, initMNS.mean) and
% quantile(data, 1-initMNS.mean)
% 'fix' - means as defined in initMNS.mean;
initMNS.mean = [];
initMNS.stdType = 'fix';
%options: 'fix', 'one group', 'multiple group', 'one group - factor', 'multiple group - factor'
%'one group - data - factor'
initMNS.temporalClasses = 'global'; % are the obsPotentialFunction defined for all time frames similar (global)
%or for each individual (individual) options: 'global, individual
initMNS.stdVal = 1; %Value of standard deviation if it is fixed to a group
initMNS.normObsProb = 2;
initMNS.neighProbFuncType = 1;
initMNS.initLabelsType = 'ones'; %options 'random', 'zeros', 'middle cluster', 'initLabels', 'ones 
initMNS.initLabelsTypeApplyToAll = true;%false; %If true applys to all metabolites if false only to unknowns
initMNS.biolReplicates = false;
initMNS.biolReplicatesAvg = 'median'; % options mean, median, integral (not yet implemented)
initMNS.determinePvalue = 0; % 0 no p-value determination, 1 only p-value determination on max breaking reaction, 2 p-value
initMNS.permutations = 1000;
initMNS.parallel.useCluster = true;
initMNS.parallel.clusterName = 'local';
initMNS.parallel.cores = 4;
initMNS.parallel.parallelFolder = '';
initMNS.parallel.parallelFolderCopy = '';
initMNS.verbose = false;
initMNS.verboseScan = 2;
initMNS.input.inferenceParameter = 2;
initMNS.input.inferenceAlgorithm = 'LazyFlipper';
initMNS.input.initLabels = [];
initMNS.input.initLabelsGene = [];
initMNS.input.reactionPairs = [];
%     initMNS.input.reactionPairs = [];
initMNS.input.edgeCliqueFile = '';
initMNS.input.edgeCliqueFileExists = false;
% initMNS.input.inferenceParameter = 2;
initMNS.input.dataFolder = '';
initMNS.input.mnsExecFolder = '';
initMNS.subGraphExtraction.std0Steps = [3 6 9 12]./3;
initMNS.subGraphExtraction.std1Steps = [0 1 2 3];
initMNS.subGraphExtraction.maxData = [];
initMNS.subGraphExtraction.mergeFreqCo = 0.1;
initMNS.subGraphExtraction.geneStd0Steps = [3 6 9 12]./3;
initMNS.subGraphExtraction.geneStd1Steps = [0 1 2 3];
initMNS.subGraphExtraction.geneMaxData = [];
initMNS.subGraphExtraction.geneStdVal = 1;
initMNS.subGraphExtraction.geneMean = [0 3];
initMNS.subGraphExtraction.geneAvgType = 'median'; %options:'all','mean','median','max';
initMNS.subGraphExtraction.geneObservations = 'best'; %options: 'all', 'best' (best = gene with max(abs(data))
initMNS.subGraphExtraction.noLogTransform = false;

for i = 1:2:nargin
    if isfield(initMNS, varargin{i})
        initMNS.(varargin{i}) = varargin{i+1};
    elseif isfield(initMNS.parallel, varargin{i})
        initMNS.parallel.(varargin{i}) = varargin{i+1};
    elseif isfield(initMNS.input, varargin{i})
        initMNS.input.(varargin{i}) = varargin{i+1};
    elseif isfield(initMNS.subGraphExtraction, varargin{i})
        initMNS.subGraphExtraction.(varargin{i}) = varargin{i+1};    
    elseif strcmp(varargin{i},'modeSubGraphExtraction')
        initMNS.meanType = 'fix';
        initMNS.mean = [varargin{i+1}.minData varargin{i+1}.maxData varargin{i+1}.maxData];
        initMNS.stdType = 'fix';
        
        if isfield(varargin{i+1}, 'std0')
            initMNS.subGraphExtraction.std0Steps = varargin{i+1}.std0;
        end
        
        if isfield(varargin{i+1}, 'std1')
            initMNS.subGraphExtraction.std1Steps = varargin{i+1}.std1;
        end
        
        initMNS.stdVal = [initMNS.subGraphExtraction.std0Steps(1) ...
            varargin{i+1}.maxData-initMNS.subGraphExtraction.std1Steps(1) ...
            varargin{i+1}.maxData-initMNS.subGraphExtraction.std1Steps(1)];
        initMNS.noOfclusters = 3;
        initMNS.initLabelsTypeApplyToAll = true;
        initMNS.initLabelsType = 'zeros';
        initMNS.normObsProb = 2;
        initMNS.normObsProb = 2;
        initMNS.nL1steps = 10;
        initMNS.subGraphExtraction.maxData = [varargin{i+1}.maxData];
    elseif strcmp(varargin{i},'modeSubGraphExtractionMetGene')
        initMNS.meanType = 'fix';
        initMNS.mean = [varargin{i+1}.minData varargin{i+1}.maxData varargin{i+1}.maxData];
        initMNS.subGraphExtraction.geneMean = [varargin{i+1}.minGeneData varargin{i+1}.maxGeneData varargin{i+1}.maxGeneData];
        initMNS.stdType = 'fix';
        
        if isfield(varargin{i+1}, 'std0')
            initMNS.subGraphExtraction.std0Steps = varargin{i+1}.std0;
        end
        
        if isfield(varargin{i+1}, 'std1')
            initMNS.subGraphExtraction.std1Steps = varargin{i+1}.std1;
        end
        if isfield(varargin{i+1}, 'geneStd0')
            initMNS.subGraphExtraction.geneStd0Steps = varargin{i+1}.geneStd0;
        end
        if isfield(varargin{i+1}, 'geneStd1')
            initMNS.subGraphExtraction.geneStd1Steps = varargin{i+1}.geneStd1;
        end
        
        initMNS.stdVal = [initMNS.subGraphExtraction.std0Steps(1) ...
            varargin{i+1}.maxData-initMNS.subGraphExtraction.std1Steps(1) ...
            varargin{i+1}.maxData-initMNS.subGraphExtraction.std1Steps(1)];
        initMNS.noOfclusters = 3;
        initMNS.initLabelsTypeApplyToAll = true;
        initMNS.initLabelsType = 'zeros';
        initMNS.normObsProb = 2;
        initMNS.nL1steps = 20;
        initMNS.subGraphExtraction.maxData = [varargin{i+1}.maxData];
        initMNS.subGraphExtraction.geneMaxData = [varargin{i+1}.maxGeneData];
        initMNS.subGraphExtraction.geneStdVal = [initMNS.subGraphExtraction.geneStd0Steps(1) ...
            varargin{i+1}.maxGeneData-initMNS.subGraphExtraction.geneStd1Steps(1) ...
            varargin{i+1}.maxGeneData-initMNS.subGraphExtraction.geneStd1Steps(1)];
    elseif strcmp(varargin{i},'temporalModel')
        initMNS.input.inferenceParameter = 1;
        initMNS.nL1steps = 20;
        initMNS.tL3steps = 20;
%         initMNS.stdType = 'one group - all data';
    end
    
end


end