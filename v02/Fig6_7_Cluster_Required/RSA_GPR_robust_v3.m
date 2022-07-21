% Compare likelihood between model structures for a variety of parameter
% sets to assess robustness of each model Models:
% 1. RSA gate and local reflex
% 2. RSA gate only
% 3. Local reflex only
% 4. No local reflex, no RSA gate
% Analysis to compare models using gaussian process regression where the
% respiratory signal is the input variable and the model is the output
% variable. The model with the highest variance is the best RSA model
% because more variance means that a change in the input is causing a
% change in the output

clear; close all;
%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);

%% Create parameter set and number of simulations
Params916030 = [2.798268 49.889025 20.086178 6.873361 0.993574 2.431832 284.973370 3.763416 1.835671 7.793137 14.090266 3.085521 14.572550 10.458438 0.183496];
% Params918720 = [2.419183 15.060474 18.687518 8.725704 0.599846 16.389965 436.648341 2.662578 2.249942 5.569219 12.253490 7.706458 6.761097 6.653920 0.973047]; %918720
% Params918723 = [2.684739 3.909452 18.148635 9.669382 3.858773 2.475236 221.908085 2.966627 1.151354 5.727808 8.894146 3.170234 3.286554 15.885458 0.564514];
n = 100; % number of sobol sample sets]
% LowerBound = [1.8192 1.5605 11.1875 1.2257 0.0998 1.7900 226.6483 0.4626 1.9999 2.5692 8.2535 2.2065 0.7611 1.1539 0.9];
LowerBound916030 = [1.8983
0.2983
15.0862
3.8734
0.1936
0.4318
34.9734
0.7634
0.6357
1.2931
2.0903
0.5855
1.45
1.4584
0.0335]'; 
Difference = Params916030 - LowerBound916030;
UpperBound916030 = Params916030 + Difference;


[sobolSet] = SobolSample(LowerBound916030, UpperBound916030,15,n);
sampleSet = sobolSet;
numSims = n;

%% RSA gate and local reflex
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13';

% Create an array of simulation input objects and specify the sweep value for each simulation

simIn(1:numSims) = Simulink.SimulationInput(mdlName);
for idx = 1:numSims
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);
    
    simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(sampleSet(idx,1)), ...        % Set ICN parameters
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(sampleSet(idx,2)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(sampleSet(idx,3)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(sampleSet(idx,4)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(sampleSet(idx,5)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(sampleSet(idx,6)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(sampleSet(idx,7)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(sampleSet(idx,8)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(sampleSet(idx,9)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(sampleSet(idx,10)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(sampleSet(idx,11)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(sampleSet(idx,12)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(sampleSet(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(sampleSet(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,15)));
end

%% Run simulations
simOut = parsim(simIn);
% simOut = sim(simIn(1,3));

%% Fit model, calculate residuals, variance, loss
filename = 'RSA_GPR_base_100_HR10916030_v13.mat';

% Analysis
%start & end times for calculating steady state values
tmeasure = [194.8 200];
tplot = [150 180];
ttest = [135 149.99];

L = zeros(1,n);
% resid = zeros(1,n);
var_base = zeros(1,n);

% numSims = 1;
j = 1;
for i = 1:numSims
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
  
    if isempty(tplotIdx) 
        IdxNAN(j) = i;
        L(1,i) = nan;
        var_base(1,i) = nan;
        j = j + 1;
    end
end

simIdx = 1:1:100;
IdxVec = setxor(simIdx,IdxNAN);

for i = 1:length(IdxVec)
    simTime         = simOut(1,IdxVec(i)).time;
    size(simTime)
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
    tmeasureIdx     = find(simTime>= tmeasure(1)  & simTime <= tmeasure(2));
    ttestIdx        = find(simTime>= ttest(1)  & simTime <= ttest(2));
    X = simOut(1,IdxVec(i)).Vlung(tplotIdx); % predictor, lung volume
    yboth = simOut(1,IdxVec(i)).HR(tplotIdx);
    gprMdl = fitrgp(X,yboth,'Standardize',true,'BasisFunction','linear');
    % Default kernel function worked best,'KernelFunction','rationalquadratic');
    % Regularization didn't decrease loss,'Regularization',0.5);
    % linear basis function has lower loss compared to quadratic and constant
    
    % make predictions and calculate loss
    ypred = resubPredict(gprMdl);
    xtest = simOut(1,IdxVec(i)).Vlung(ttestIdx);
    ytest = simOut(1,IdxVec(i)).HR(ttestIdx);
    L(1,IdxVec(i)) = loss(gprMdl,xtest,ytest);
    disp(IdxVec(i))
    res = simOut(1,IdxVec(i)).HR(tplotIdx) - ypred;
    mean(res)
    length(res)
    var_base(1,IdxVec(i)) = var(res);
    
end
save(filename,'L','var_base')

%% RSA gate only
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13_NLR';

simIn(1:numSims) = Simulink.SimulationInput(mdlName);
for idx = 1:numSims
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);
    
    simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(sampleSet(idx,1)), ...        % Set ICN parameters
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(sampleSet(idx,2)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(sampleSet(idx,3)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(sampleSet(idx,4)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(sampleSet(idx,5)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(sampleSet(idx,6)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(sampleSet(idx,7)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(sampleSet(idx,8)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(sampleSet(idx,9)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(sampleSet(idx,10)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(sampleSet(idx,11)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(sampleSet(idx,12)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(sampleSet(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(sampleSet(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,15)));
end

% Run simulations
simOut = parsim(simIn);

% Fit model, calculate residuals, variance, loss
filename = 'RSA_GPR_NLR_100_HR10916030_v13.mat';

% Analysis
%start & end times for calculating steady state values
tmeasure = [194.8 200];
tplot = [150 179.9];
ttest = [135 149.99];

L = zeros(1,n);
% resid = zeros(1,n);
var_base = zeros(1,n);

% numSims = 1;
j = 1;
for i = 1:numSims
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
  
    if isempty(tplotIdx) 
        IdxNAN(j) = i;
        L(1,i) = nan;
        var_base(1,i) = nan;
        j = j + 1;
    end
end

simIdx = 1:1:100;
IdxVec = setxor(simIdx,IdxNAN);

for i = 1:length(IdxVec)
    simTime         = simOut(1,IdxVec(i)).time;
    size(simTime)
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
    tmeasureIdx     = find(simTime>= tmeasure(1)  & simTime <= tmeasure(2));
    ttestIdx        = find(simTime>= ttest(1)  & simTime <= ttest(2));
    X = simOut(1,IdxVec(i)).Vlung(tplotIdx); % predictor, lung volume
    yboth = simOut(1,IdxVec(i)).HR(tplotIdx);
    gprMdl = fitrgp(X,yboth,'Standardize',true,'BasisFunction','linear');
    % Default kernel function worked best,'KernelFunction','rationalquadratic');
    % Regularization didn't decrease loss,'Regularization',0.5);
    % linear basis function has lower loss compared to quadratic and constant
    
    % make predictions and calculate loss
    ypred = resubPredict(gprMdl);
    xtest = simOut(1,IdxVec(i)).Vlung(ttestIdx);
    ytest = simOut(1,IdxVec(i)).HR(ttestIdx);
    L(1,IdxVec(i)) = loss(gprMdl,xtest,ytest);
    disp(IdxVec(i))
    res = simOut(1,IdxVec(i)).HR(tplotIdx) - ypred;
    mean(res)
    length(res)
    var_base(1,IdxVec(i)) = var(res);
    
end
save(filename,'L','var_base')

%% Local reflex only
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13_noRSA';

simIn(1:numSims) = Simulink.SimulationInput(mdlName);
for idx = 1:numSims
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);
    
    simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(sampleSet(idx,1)), ...        % Set ICN parameters
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(sampleSet(idx,2)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(sampleSet(idx,3)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(sampleSet(idx,4)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(sampleSet(idx,5)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(sampleSet(idx,6)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(sampleSet(idx,7)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(sampleSet(idx,8)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(sampleSet(idx,9)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(sampleSet(idx,10)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(sampleSet(idx,11)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(sampleSet(idx,12)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(sampleSet(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(sampleSet(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,15)));
end

% Run simulations
simOut = parsim(simIn);

% Fit model, calculate residuals, variance, loss
filename = 'RSA_GPR_noRSA_100_HR10916030_v13.mat';


% Analysis
%start & end times for calculating steady state values
tmeasure = [194.8 200];
tplot = [150 179.9];
ttest = [135 149.99];

L = zeros(1,n);
% resid = zeros(1,n);
var_base = zeros(1,n);

% numSims = 1;
j = 1;
for i = 1:numSims
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
  
    if isempty(tplotIdx) 
        IdxNAN(j) = i;
        L(1,i) = nan;
        var_base(1,i) = nan;
        j = j + 1;
    end
end

simIdx = 1:1:100;
IdxVec = setxor(simIdx,IdxNAN);

for i = 1:length(IdxVec)
    simTime         = simOut(1,IdxVec(i)).time;
    size(simTime)
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
    tmeasureIdx     = find(simTime>= tmeasure(1)  & simTime <= tmeasure(2));
    ttestIdx        = find(simTime>= ttest(1)  & simTime <= ttest(2));
    X = simOut(1,IdxVec(i)).Vlung(tplotIdx); % predictor, lung volume
    yboth = simOut(1,IdxVec(i)).HR(tplotIdx);
    gprMdl = fitrgp(X,yboth,'Standardize',true,'BasisFunction','linear');
    % Default kernel function worked best,'KernelFunction','rationalquadratic');
    % Regularization didn't decrease loss,'Regularization',0.5);
    % linear basis function has lower loss compared to quadratic and constant
    
    % make predictions and calculate loss
    ypred = resubPredict(gprMdl);
    xtest = simOut(1,IdxVec(i)).Vlung(ttestIdx);
    ytest = simOut(1,IdxVec(i)).HR(ttestIdx);
    L(1,IdxVec(i)) = loss(gprMdl,xtest,ytest);
    disp(IdxVec(i))
    res = simOut(1,IdxVec(i)).HR(tplotIdx) - ypred;
    mean(res)
    length(res)
    var_base(1,IdxVec(i)) = var(res);
    
end
save(filename,'L','var_base')

%% No local reflex no RSA
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13_NLR_NRSA';

simIn(1:numSims) = Simulink.SimulationInput(mdlName);
for idx = 1:numSims
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);
    
    simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(sampleSet(idx,1)), ...        % Set ICN parameters
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(sampleSet(idx,2)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(sampleSet(idx,3)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(sampleSet(idx,4)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(sampleSet(idx,5)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(sampleSet(idx,6)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(sampleSet(idx,7)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(sampleSet(idx,8)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(sampleSet(idx,9)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(sampleSet(idx,10)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(sampleSet(idx,11)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(sampleSet(idx,12)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(sampleSet(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(sampleSet(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,15)));
end

% Run simulations
simOut = parsim(simIn);

% Fit model, calculate residuals, variance, loss
filename = 'RSA_GPR_NLR_NRSA_100_HR10916030_v13.mat';


% Analysis
%start & end times for calculating steady state values
tmeasure = [194.8 200];
tplot = [150 179.9];
ttest = [135 149.99];

L = zeros(1,n);
% resid = zeros(1,n);
var_base = zeros(1,n);

% numSims = 1;
j = 1;
for i = 1:numSims
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
  
    if isempty(tplotIdx) 
        IdxNAN(j) = i;
        L(1,i) = nan;
        var_base(1,i) = nan;
        j = j + 1;
    end
end

simIdx = 1:1:100;
IdxVec = setxor(simIdx,IdxNAN);

for i = 1:length(IdxVec)
    simTime         = simOut(1,IdxVec(i)).time;
    size(simTime)
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
    tmeasureIdx     = find(simTime>= tmeasure(1)  & simTime <= tmeasure(2));
    ttestIdx        = find(simTime>= ttest(1)  & simTime <= ttest(2));
    X = simOut(1,IdxVec(i)).Vlung(tplotIdx); % predictor, lung volume
    yboth = simOut(1,IdxVec(i)).HR(tplotIdx);
    gprMdl = fitrgp(X,yboth,'Standardize',true,'BasisFunction','linear');
    % Default kernel function worked best,'KernelFunction','rationalquadratic');
    % Regularization didn't decrease loss,'Regularization',0.5);
    % linear basis function has lower loss compared to quadratic and constant
    
    % make predictions and calculate loss
    ypred = resubPredict(gprMdl);
    xtest = simOut(1,IdxVec(i)).Vlung(ttestIdx);
    ytest = simOut(1,IdxVec(i)).HR(ttestIdx);
    L(1,IdxVec(i)) = loss(gprMdl,xtest,ytest);
    disp(IdxVec(i))
    res = simOut(1,IdxVec(i)).HR(tplotIdx) - ypred;
    mean(res)
    length(res)
    var_base(1,IdxVec(i)) = var(res);
    
end
save(filename,'L','var_base')