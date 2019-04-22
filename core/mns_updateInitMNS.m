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
function initMNS = mns_updateInitMNS(initMNS, varargin)

if nargin < 1 || isempty(initMNS)
    initMNS = mns_generateInitMNS;
end

for i = 1:2:nargin-1
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
        initMNS.nL1steps = 20;
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
        initMNS.normObsProb = 2;
        initMNS.nL1steps = 10;
        initMNS.subGraphExtraction.maxData = [varargin{i+1}.maxData];
        initMNS.subGraphExtraction.geneMaxData = [varargin{i+1}.maxGeneData];
        initMNS.subGraphExtraction.geneStdVal = [initMNS.subGraphExtraction.geneStd0Steps(1) ...
            varargin{i+1}.maxGeneData-initMNS.subGraphExtraction.geneStd1Steps(1) ...
            varargin{i+1}.maxGeneData-initMNS.subGraphExtraction.geneStd1Steps(1)];
    end
    
end

end