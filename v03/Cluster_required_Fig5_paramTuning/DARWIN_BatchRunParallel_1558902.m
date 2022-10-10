% Script to tune ICN parameters using data from Rajendran 2019
% Heart rate and firing frequency changes calculated from Rajendran 2019
% Calculations can be seen in Fig7_VNS.m
% Data set from Pennsieve: https://app.pennsieve.io/N:organization:618e8dd9-f8d2-4dc4-9abb-c6aaab2e78a0/datasets/N:dataset:0d9454ca-43d9-4fdb-ac79-01c3574e8565/files/N:collection:063f029d-ebf1-49fc-803c-1c618f314ba4
% Protocol: https://www.protocols.io/view/pig-icn-recording-2jugcnw?step=4

clear;
%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);


%% Define name of simulink model
mdlName     = 'ICN_model_v15_tune';

%% Create parameter set and number of simulations
n = 5000; % number of sobol sample sets
[sobolSet] = SobolSample([0.06,1.3,2,0.1,0.06,2.3,2,0.1,0.06,5,5,11,1,1,1,0.001], [8,75,22,8,10,150,700,5,3,20,17,14,10,10,10,0.999],16,n);
sampleSet = sobolSet;
numSims = n;
PNDMVdelay = 0.3;

%% Create an array of simulation input objects and specify the sweep value for each simulation

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
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(PNDMVdelay), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(sampleSet(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(sampleSet(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(sampleSet(idx,15)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,16)));
end

%% Run simulations
simout = parsim(simIn);
% simout(1,1) = sim(simIn(1,1));
%% Analysis
avgHRpre= zeros(1,n);

x = [1.59;1.54;1.62;1.59;1.54;1.62;1.59;1.54;1.62;67.5503;2.887;-0.0767*880;0.0018*37528];% tuning targets using poisson dist to estimate
% old tuning targets [1.59;1.57;1.62;1.59;1.57;1.62;1.59;1.57;1.62;67.5503;2.887];%84.6;2.392]; % Baseline, stimulation, recovery period averaged ICN neuron response to vagal stimulation, hR and Emaxlv from Ursino 1998, which park used for fitting
% Average across 7 neurons recorded in Rajendran data was used because all
% ICN neurons should roughly fire at this frequency
% Numbers are repeated three times, for NAPN,DMVPN, and LCN
% Park model output HR and Emaxlv: 67.5503;2.8871
% See Fig7_VNS.m for analysis of Rajendran 2019 data to obtain
% percentrage changes in HR in response to 1 Hz step input

% HRpercentchangeduringstimulation =
% 
%    -0.0767
% 
% 
% HRpercentchangeduring recovery =
% 
%     0.0018
    
simRes = zeros(length(x),n);
mse = zeros(1,n);

for i = 1:numSims
    
    failedSim = isnan(avgHRpre);
    if failedSim == 1
        for j = 1:length(x)
            simRes(j,i) = NaN;
        end
    else

        HR_vals         = simout(1,i).HR_TS;
        Emaxlv_vals     = simout(1,i).Emaxlv;

        
        PN_NA_ff	= simout(1,i).PN_NA_ff;
        PN_DMV_ff	= simout(1,i).PN_DMV_ff;
        LCN_ff      = simout(1,i).LCN_ff;
        
        % Calculate percent change in heart rate due to step
        if length(mean(getsampleusingtime(HR_vals,240,300),'Weighting','time')) == 0
            simRes(1,i) 	= nan; %baseline
            simRes(2,i)     = nan; %VNS
            simRes(3,i) 	= nan; % post VNS
            simRes(4,i)     = nan; %baseline
            simRes(5,i)     = nan; %VNS
            simRes(6,i)     = nan; % post VNS
            simRes(7,i)     = nan; %baseline
            simRes(8,i)     = nan; %VNS
            simRes(9,i)     = nan; % post VNS
            simRes(10,i)    = nan;
            simRes(11,i)    = nan;
            simRes(12,i)    = nan;
            simRes(13,i)    = nan;
        else
            initialHR = mean(getsampleusingtime(HR_vals,165,180),'Weighting','time');
            initialEmax = mean(getsampleusingtime(Emaxlv_vals,165,180),'Weighting','time');
            simRes(1,i) 	= mean(getsampleusingtime(PN_NA_ff,120,180),'Weighting','time'); %baseline
            simRes(2,i)     = mean(getsampleusingtime(PN_NA_ff,180,240),'Weighting','time'); %VNS
            simRes(3,i) 	= mean(getsampleusingtime(PN_NA_ff,240,300),'Weighting','time'); % post VNS
            simRes(4,i)     = mean(getsampleusingtime(PN_DMV_ff,120, 180),'Weighting','time'); %baseline
            simRes(5,i)     = mean(getsampleusingtime(PN_DMV_ff,180,240),'Weighting','time'); %VNS
            simRes(6,i)     = mean(getsampleusingtime(PN_DMV_ff,240,300),'Weighting','time'); % post VNS
            simRes(7,i)     = mean(getsampleusingtime(LCN_ff,120, 180),'Weighting','time'); %baseline
            simRes(8,i)     = mean(getsampleusingtime(LCN_ff,180,240),'Weighting','time'); %VNS
            simRes(9,i)     = mean(getsampleusingtime(LCN_ff,240,300),'Weighting','time'); % post VNS
            simRes(10,i)    = initialHR;
            simRes(11,i)    = initialEmax;
            simRes(12,i)    = ((mean(getsampleusingtime(HR_vals,180,240),'Weighting','time') - initialHR) / initialHR) *880;
            simRes(13,i)    = ((mean(getsampleusingtime(HR_vals,240,300),'Weighting','time') - initialHR) / initialHR) *37528;
        end
    end
    
end


for i = 1:numSims
    noVal = isnan(simRes(1,i));
    
    %calculate mse if all simRes components were calculated
    if noVal == 1
        mse(i) = 100000;
    else
        mse(i) = immse(x,simRes(:,i));
    end
end


%% find 10 parameter sets with lowest mse
minMSEparamSets = 10;
paramMatrix = zeros((minMSEparamSets),length(sampleSet(1,:)));
minidxstore = zeros(1,minMSEparamSets);
for i = 1:minMSEparamSets
    [msemin, minidx]  = min(mse)
    minidxstore(i) = minidx;
    disp('Parameters used: PN-NA {fmin, fmax, fmid, gain}, LCN {fmin, fmax, fmid, gain}, PN-DMV {fmin, fmax, fmid, gain}, kfev, kBR, kCP, kfesh')
    sprintf('% 6.6f',sampleSet(minidx,:)')
    disp('Initial HR')
    errormessage = simout(1,i).ErrorMessage
    mse(minidx) = 100000.1001;
    paramMatrix(i,:) = sampleSet(minidx,:);
end


%% Calculate min and max of each parameter to determine range for next sobol
% sampling round
minParam = zeros(1,length(sampleSet(1,:)));
maxParam = zeros(1,length(sampleSet(1,:)));
for i = 1:length(sampleSet(1,:))
    minParam(i) = min(paramMatrix(:,i));
    maxParam(i) = max(paramMatrix(:,i));
end
sprintf('% 6.6f',minParam')
sprintf('% 6.6f',maxParam')


%% Second round of sobol sampling
disp('Second round of sobol sampling')

% Create parameter set and number of simulations
[sobolSet] = SobolSample(minParam', maxParam', 16, n);
sampleSet = sobolSet;
numSims = n;

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
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(PNDMVdelay), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(sampleSet(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(sampleSet(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(sampleSet(idx,15)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,16)));
end

%% Run simulations
simout = parsim(simIn);

%% Analysis
% Heart rates during different time periods
avgHRpre= zeros(1,n);

x = [1.59;1.54;1.62;1.59;1.54;1.62;1.59;1.54;1.62;67.5503;2.887;-0.0767;0.0018]; % tuning targets using poisson dist to estimate
% old tuning targets [1.59;1.57;1.62;1.59;1.57;1.62;1.59;1.57;1.62;67.5503;2.887];%84.6;2.392]; % Baseline, stimulation, recovery period averaged ICN neuron response to vagal stimulation, hR and Emaxlv from Ursino 1998, which park used for fitting
% Average across 7 neurons recorded in Rajendran data was used because all
% ICN neurons should roughly fire at this frequency
% Numbers are repeated three times, for NAPN,DMVPN, and LCN
% Park model output HR and Emaxlv: 67.5503;2.8871
simRes = zeros(length(x),n);
mse = zeros(1,n);

for i = 1:numSims

    failedSim = isnan(avgHRpre);
    if failedSim == 1
        for j = 1:length(x)
            simRes(j,i) = NaN;
        end
    else

        HR_vals         = simout(1,i).HR_TS;
        Emaxlv_vals     = simout(1,i).Emaxlv;
        
        PN_NA_ff	= simout(1,i).PN_NA_ff;
        PN_DMV_ff	= simout(1,i).PN_DMV_ff;
        LCN_ff      = simout(1,i).LCN_ff;

        
        % Calculate percent change in heart rate due to step
        if length(mean(getsampleusingtime(HR_vals,240,300),'Weighting','time')) == 0
            simRes(1,i) 	= nan; %baseline
            simRes(2,i)     = nan; %VNS
            simRes(3,i) 	= nan; % post VNS
            simRes(4,i)     = nan; %baseline
            simRes(5,i)     = nan; %VNS
            simRes(6,i)     = nan; % post VNS
            simRes(7,i)     = nan; %baseline
            simRes(8,i)     = nan; %VNS
            simRes(9,i)     = nan; % post VNS
            simRes(10,i)    = nan;
            simRes(11,i)    = nan;
            simRes(12,i)    = nan;
            simRes(13,i)    = nan;
        else
            initialHR = mean(getsampleusingtime(HR_vals,165,180),'Weighting','time');
            initialEmax = mean(getsampleusingtime(Emaxlv_vals,165,180),'Weighting','time');
            simRes(1,i) 	= mean(getsampleusingtime(PN_NA_ff,120,180),'Weighting','time'); %baseline
            simRes(2,i)     = mean(getsampleusingtime(PN_NA_ff,180,240),'Weighting','time'); %VNS
            simRes(3,i) 	= mean(getsampleusingtime(PN_NA_ff,240,300),'Weighting','time'); % post VNS
            simRes(4,i)     = mean(getsampleusingtime(PN_DMV_ff,120, 180),'Weighting','time'); %baseline
            simRes(5,i)     = mean(getsampleusingtime(PN_DMV_ff,180,240),'Weighting','time'); %VNS
            simRes(6,i)     = mean(getsampleusingtime(PN_DMV_ff,240,300),'Weighting','time'); % post VNS
            simRes(7,i)     = mean(getsampleusingtime(LCN_ff,120, 180),'Weighting','time'); %baseline
            simRes(8,i)     = mean(getsampleusingtime(LCN_ff,180,240),'Weighting','time'); %VNS
            simRes(9,i)     = mean(getsampleusingtime(LCN_ff,240,300),'Weighting','time'); % post VNS
            simRes(10,i)    = initialHR;
            simRes(11,i)    = initialEmax;
            simRes(12,i)    = (mean(getsampleusingtime(HR_vals,180,240),'Weighting','time') - initialHR) / initialHR;
            simRes(13,i)    = (mean(getsampleusingtime(HR_vals,240,300),'Weighting','time') - initialHR) / initialHR;
        end
    end

end

for i = 1:numSims
    noVal = isnan(simRes(1,i));
    
    %calculate mse if all simRes components were calculated
    if noVal == 1
        mse(i) = 100000;
    else
        mse(i) = immse(x,simRes(:,i));
    end
end


%% find 10 parameter sets with lowest mse
minMSEparamSets = 10;
paramMatrix = zeros((minMSEparamSets),length(sampleSet(1,:)));
minidxstore = zeros(1,minMSEparamSets);
for i = 1:minMSEparamSets
    [msemin, minidx]  = min(mse)
    minidxstore(i) = minidx;
    disp('Parameters used: PN-NA {fmin, fmax, fmid, gain}, LCN {fmin, fmax, fmid, gain}, PN-DMV {fmin, fmax, fmid, gain}, kfev, kBR, kCP, kfesh')
    sprintf('% 6.6f',sampleSet(minidx,:)')
    disp('Initial HR')
    errormessage = simout(1,i).ErrorMessage
    mse(minidx) = 100000.1001;
    paramMatrix(i,:) = sampleSet(minidx,:);
end

%% Exit code
close_system
delete(myPool);
exit