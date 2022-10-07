% Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
% Michelle Gee
% October 3, 2022

% Script to produce Figure 5

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

%% Create simulation inputs
Params = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 3.329861 2.661685 5.642977 0.066794];
LowerBound = 0;
UpperBound = 1;
kRSA = LowerBound:0.1:UpperBound;

numSims = length(kRSA);
n = numSims;

mdlName = 'ICN_model_v15_VNS_cluster';

% Create an array of simulation input objects and specify the sweep value for each simulation
simIn(1:numSims) = Simulink.SimulationInput(mdlName);
for idx = 1:numSims
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);
    
    simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...        % Set ICN parameters
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(Params(2)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(Params(3)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(Params(4)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(Params(5)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(Params(6)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(Params(7)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(Params(8)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(Params(9)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(Params(10)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(Params(11)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(Params(12)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(Params(15)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(kRSA(idx)));
end

%% Run simulations
simOut = parsim(simIn);
% simOut = sim(simIn(1,1));

%% Fit model, calculate residuals, variance, loss
filename = 'RSA_GPR_kRSA_09_30_22.mat';

% Analysis
%start & end times for calculating steady state values
tplot = [150 179.9];
ttest = [135 149.99];

% preallocate loss and variance vectors for storage
L = zeros(1,n);
var_base = zeros(1,n);

% determine indices for each simulation
j = 1;
for i = 1:numSims
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
  
    % check if simulation failed
    if isempty(tplotIdx) 
        IdxNAN(j) = i;
        L(1,i) = nan;
        var_base(1,i) = nan;
        j = j + 1;
    end
end

% define simulation number
simIdx = 1:1:length(kRSA);
IdxVec = simIdx;

for i = 1:length(IdxVec)
    simTime         = simOut(1,IdxVec(i)).time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2)); % training indices
    ttestIdx        = find(simTime>= ttest(1)  & simTime <= ttest(2)); % test indices
    X = simOut(1,IdxVec(i)).Vlung(tplotIdx)'; % predictor, lung volume, must be a column vector for normalize
    X_norm = normalize(X);  % normalize data
    X_row = X_norm'; % row vector for fitrgp
    yboth = simOut(1,IdxVec(i)).HR(tplotIdx)'; % HR, value to be predicted
    y_norm = normalize(yboth); % normalize data
    y_row = y_norm';
    gprMdl = fitrgp(X_row,y_row,'BasisFunction','linear');
    % Default kernel function worked best,'KernelFunction','rationalquadratic');
    % Regularization didn't decrease loss,'Regularization',0.5);
    % linear basis function has lower loss compared to quadratic and constant
    
    % make predictions
    ypred = resubPredict(gprMdl); % predict HR from model
    xtest = simOut(1,IdxVec(i)).Vlung(ttestIdx)'; % extract and normalize test data
    xtest_norm = normalize(xtest)';
    ytest = simOut(1,IdxVec(i)).HR(ttestIdx)';
    ytest_norm = normalize(ytest)';
    
    % calculate loss
    L(1,IdxVec(i)) = loss(gprMdl,xtest_norm,ytest_norm);
    mdlHR_norm = normalize(simOut(1,IdxVec(i)).HR(tplotIdx)')'; % normalize along cols, but gpr regression produces row vectors
    
    % calculate residual and variance
    res = mdlHR_norm - ypred;
    var_base(1,IdxVec(i)) = var(res);
    
end

% save residuals and variances for each model
save(filename,'L','var_base')
% plot(simTime(tplotIdx),ypred,'o',simTime(tplotIdx),mdlHR_norm,'--')
