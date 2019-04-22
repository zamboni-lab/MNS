%% Example Identification of sites and sequential order from metabolomics data via metabolic network segmentation

%% Summary:
% example code to infer sites and sequential order of metabolic regulation
% in fibroblasts treated with different concentrations of H2O2 (Data from
% Kuehne et al)

%% load the data
load('mns_example_sequential_mns_fibroblasts_H2O2 - WS.mat')


%% initialize mns
% go to mns folder and perform the execute command:
mns_initialize

%% perform automated coarse grained MNS inference
% initialize mns parameter settings
initMNS = mns_generateInitMNS('inferenceParameter',1,'nL1steps',20, 'tL3steps', 20);
initMNS.noOfclusters = 3; 
initMNS.mean = [-0.10 0 0.10];
initMNS.meanType = 'fix';
initMNS.stdType = 'one group - all data';

model = KEGG_HSA_MNS_red;

% run the mns algorithm
mnsResultsFibroDill_coarse = mns_scanTime(model,dataStructFibroH2O2, 'fibro_dil_time_coarse',0,0,initMNS);
%% plot the fracture frequency and score distribution

% fracture frequency and sum of observation potential
mns_plotFractureFrequency(mnsResultsFibroDill_coarse, model,0.25,0.25,0.01,0.01)

% Score distributions for increasing weights of ws and wt
mns_plotMultipleScreeningResults(mnsResultsFibroDill_coarse, model,0.25,0.25,0.05,0.05)

% Module label, sequence fracture and neighborhodo fracture distribution
% for wn = ws = 0.0 (no influence)
mns_calcProbabilityScanData(mnsResultsFibroDill_coarse, model,0.0,0.0, false);

% Module label, sequence fracture and neighborhodo fracture distribution
% for wn = ws = 0.2
[~,~,idx] = mns_calcProbabilityScanData(mnsResultsFibroDill_coarse, model,0.21,0.21, false);

% Module label, sequence fracture and neighborhodo fracture distribution
% for wn = ws = 0.25
mns_calcProbabilityScanData(mnsResultsFibroDill_coarse, model,0.25,0.25, false);
colormap([0 102 204; 255 255 255; 204 0 51]/255)
% Comment: For ws and wt > 0.21 only one sequential and one neighborhood
% fracture is left which is biologically not relevant. Therefore we used
% as ws = wn = 0 as reference point and investigated the fracture stability
% in a more fine grained scan thorugh lambda1 and lambda2 in the range from
% 0 to approximately double the maximal lambda value with the ws = wn = 0.

% Lambda values for ws = wn = 0.21;
disp(['lambda1 = ' num2str(mnsResultsFibroDill_coarse.nL1(idx))]);
disp(['lambda2 = ' num2str(mnsResultsFibroDill_coarse.tL3(idx))]);

%% perform the MNS inference
% set the scanning range of lambda1 (nlRangeTemp) and lambda2 (tlRangeTemp)
% Note: the upperlimit of both lambda values have been determined manually
% before (see section before)
tlRangeTemp = [0 0.01:0.01:0.6];
nlRangeTemp = [0 0.01:0.01:0.6];

% initialize mns parameter settings
initMNS = mns_generateInitMNS('inferenceParameter',1,'nL1steps',20, 'tL3steps', 20);
initMNS.noOfclusters = 3;
initMNS.mean = [-0.10 0 0.10];
initMNS.meanType = 'fix';
initMNS.stdType = 'one group - all data';

% run the mns algorithm
mnsResultsFibroDill_fineRange = mns_scanTime(model,dataStructFibroH2O2, 'fibro_dil_time',nlRangeTemp,tlRangeTemp,initMNS);

%% Data analysis
% Sum of neighborhood fractures, sequential frame fracture and sum of observation potential
mns_calcProbabilityScanData(mnsResultsFibroDill_fineRange, model,0.00,0.00,true, false);

%% Supplementary Figure 9b
% Freqeuency of sequential frame and neighborhood fractures with increasing
% weights ws and wn
mns_plotFractureFrequency(mnsResultsFibroDill_fineRange, model,0.5,0.5,0.01,0.01)

% select the weights at which number of fractures are at pseudo
% steady state (i.e. constant over a small range of weights)
w = [0 0.05 0.1 0.15 0.24 0.3 0.4];

%% Plot score distribution for increasing weights of ws and wn as defined before.
mns_plotMultipleScreeningResults(mnsResultsFibroDill_fineRange, model,w,w)

%% Plot Module labels, sequence fractures and neighborhood fractures for all metabolites and reactions with increasing weights ws and wn. 
mns_plotModuleLabelsAndFractures(mnsResultsFibroDill_fineRange, model,w,w)

%% PLot module labels on top of sequential metabolomics data

% left: sequence weight ws = 0, neighborhood weight wn = 0
[~,~,idxFibro] =  mns_calcProbabilityScanData(mnsResultsFibroDill_fineRange, model,0,0, false, false);
mns_plotTcWithModules(mnsResultsFibroDill_fineRange, dataStructFibroH2O2, model, idxFibro, 'area')

% middle: ws = 0.15, wn = 0.15
[~,~,idxFibro] =  mns_calcProbabilityScanData(mnsResultsFibroDill_fineRange, model,0.15,0.15, false, false);
mns_plotTcWithModules(mnsResultsFibroDill_fineRange, dataStructFibroH2O2, model, idxFibro, 'area')

% right: ws = 0.24, wn = 0.24
[~,~,idxFibro] =  mns_calcProbabilityScanData(mnsResultsFibroDill_fineRange, model,0.24,0.24, false, false);
mns_plotTcWithModules(mnsResultsFibroDill_fineRange, dataStructFibroH2O2, model, idxFibro, 'area')

% not shown: ws = 0.3, wn = 0.3
[~,~,idxFibro] =  mns_calcProbabilityScanData(mnsResultsFibroDill_fineRange, model,0.3,0.3, false, false);
mns_plotTcWithModules(mnsResultsFibroDill_fineRange, dataStructFibroH2O2, model, idxFibro, 'area')

%% Generate outtables of fractures
[seqFracList, seqFracListOnlyFracPos, nFracList, nFracListSplit] ...
    = mns_scanTimeResults2table(mnsResultsFibroDill_fineRange,model, 0.025, 0.025, 0.3);