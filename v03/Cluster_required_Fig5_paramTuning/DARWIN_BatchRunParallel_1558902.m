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
% simout(1,1) = sim(simIn(1,3));
%% Analysis
%start & end times for calculating steady state values
tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
tlong  = [285 300];     % for CO and SV determination

% Fitting PN_NA to VNS
VNStimepre = [120 180]; %after steady state reached, before VNS
VNStime     = [180 240]; %during VNS
VNStimepost = [240 300];

% Heart rates during different time periods
avgHRpre= zeros(1,n);
avgHRstim = zeros(1,n);
avgHRpost = zeros(1,n);

x = [1.59;1.54;1.62;1.59;1.54;1.62;1.59;1.54;1.62;67.5503;2.887;-0.0767*880;0.0018*37528];% tuning targets using poisson dist to estimate
% old tuning targets [1.59;1.57;1.62;1.59;1.57;1.62;1.59;1.57;1.62;67.5503;2.887];%84.6;2.392]; % Baseline, stimulation, recovery period averaged ICN neuron response to vagal stimulation, hR and Emaxlv from Ursino 1998, which park used for fitting
% Average across 7 neurons recorded in Rajendran data was used because all
% ICN neurons should roughly fire at this frequency
% Numbers are repeated three times, for NAPN,DMVPN, and LCN
% Park model output HR and Emaxlv: 67.5503;2.8871
% See Fig7_VNS.m for analysis of Rajendran 2019 data to obtain
% percentrage changes in HR in response to 1 Hz step input

% HRstimDel =
% 
%    -0.0767
% 
% 
% HRrecovDel =
% 
%     0.0018
    
simRes = zeros(length(x),n);
mse = zeros(1,n);

for i = 1:numSims
    simTime         = simout(1,i).time;
    tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
    tlongIdx 	= find(simTime >= tlong(1) & simTime <= tlong(2));
    tVNSpre 	= find(simTime >= VNStimepre(1) & simTime <= VNStimepre(2));
    tVNS        = find(simTime >= VNStime(1) & simTime <= VNStime(2));
    tVNSpost 	= find(simTime >= VNStimepost(1) & simTime <= VNStimepost(2));
    
    
    failedSim = isnan(avgHRpre);
    if failedSim == 1
        for j = 1:length(x)
            simRes(j,i) = NaN;
        end
    else
        simTime         = simout(1,i).time;
        tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
        tlongIdx        = find(simTime>= tlong(1)  & simTime <= tlong(2));
        
        Tperiod_vals    = simout(1,i).Tperiod;
        HR_vals         = simout(1,i).HR_TS;
        Emaxlv_vals     = simout(1,i).Emaxlv;
        
        Rsp_vals        = simout(1,i).Rsp;
        Rep_vals        = simout(1,i).Rep;
        Rmp_vals        = simout(1,i).Rmp;
        
        Vlv_vals        = simout(1,i).Vlv;
        
        Psa_vals        = simout(1,i).Psa;
        
        Flow_lv_vals    = simout(1,i).Flow_lv;
        Flow_rv_vals    = simout(1,i).Flow_rv;
        
        Pla_vals        = simout(1,i).Pla;
        Plv_vals        = simout(1,i).Plv;
        
        phi_vals        = simout(1,i).Phi;
        
        PN_NA_ff	= simout(1,i).PN_NA_ff;
        PN_DMV_ff	= simout(1,i).PN_DMV_ff;
        LCN_ff      = simout(1,i).LCN_ff;
        fevHR		= simout(1,i).fevHR;
        fevEmax		= simout(1,i).fevEmax;
        
        try
            CO_SV = fx_CO_SV(tlong, simTime, Flow_lv_vals, Flow_rv_vals, Tperiod_vals);
        catch
            %warning('error in CO/SV calculation likely due to parameters used');
            CO_SV = struct('RCO', NaN, 'RSV', NaN, 'LCO', NaN,  'LSV', NaN);
        end
        
        EDVthreshold    = 10; % end diastolic volume min threshold
        try
            [Vlv_pks, Vlv_index]    = findpeaks(Vlv_vals(tshortIdx), 'MinPeakHeight', EDVthreshold,  'MinPeakProminence', 5);   % returns peak values and indices of peaks
        catch
            warning('error in left ventricular volume calculation')
            Vlv_pks     = 999;
        end
        
        meanVlvPks      = mean(Vlv_pks);
        SysThreshold    = 10;
        
        try
            [Psys_pks, Psys_index]  = findpeaks(Psa_vals(tlongIdx),    'MinPeakHeight', SysThreshold, 'MinPeakProminence', 5);
            [Pdia_pks, Pdia_index]  = findpeaks(-1*Psa_vals(tlongIdx), 'MinPeakProminence', 5);
            
            if length(Psys_pks) ~= length(Pdia_pks)
                minLength   = min([length(Psys_pks), length(Pdia_pks)]);
                Psys_pks    = Psys_pks(1:minLength);
                Pdia_pks    = Pdia_pks(1:minLength);
            end
        catch
            Psys_pks    = 999*ones(10, 1);  % arbitrarily chosen vector length of 10
            Pdia_pks    = 999*ones(10, 1);  % arbitrarily chosen vector length of 10
        end
        
        % identify appropriate indices for Left Ventricualr End Systolic Volume
        % (ESV) and End Systolic Pressure (ESP)
        Vlv_vals_tshort     = Vlv_vals(tshortIdx);
        phi_vals_tshort     = phi_vals(tshortIdx);
        Vlv_idx             = find(Vlv_vals_tshort <= (meanVlvPks-10));
        phi_idx             = find(phi_vals_tshort <= 0.2);
        int_idx             = intersect(Vlv_idx, phi_idx);
        ESV_idx             = int_idx(find(diff(int_idx)> 5) + 1);
        ESV_vals            = Vlv_vals_tshort(ESV_idx);
        
        Plv_vals_tshort     = Plv_vals(tshortIdx);
        Vlv_idx2            = find(diff(Vlv_vals_tshort) <= 0.01 & diff(Vlv_vals_tshort) >= -0.01);
        int_idx2            = intersect(Vlv_idx,Vlv_idx2);
        ESP_idx             = int_idx2(find(diff(int_idx2) > 20)+1);
        ESP_vals            = Plv_vals_tshort(ESP_idx);
        
        % identify appropriate indices for Left Ventricualr End DIastolic Volume
        % (EDV) and End Diastolic Pressure (ESP)
        phi_idx             = find(phi_vals_tshort == 0);
        EDV_idx             = phi_idx(find(diff(phi_idx) > 2));
        EDV_vals            = Vlv_vals_tshort(EDV_idx);
        
        Vlv_idx3            = find(Vlv_vals_tshort > (meanVlvPks-10));
        int_idx3            = intersect(phi_idx, Vlv_idx3);
        EDP_idx             = int_idx3(find(diff(int_idx3)>2));
        EDP_vals            = Plv_vals_tshort(EDP_idx);
        
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
    disp('Parameters used: fmin, fmax, fmid, gain for NA, NActr, BR')
    sprintf('% 6.6f',sampleSet(minidx,:)')
    disp('Initial HR')
%     avgHRbase(minidx)
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

%% Create parameter set and number of simulations
[sobolSet] = SobolSample(minParam', maxParam', 16, n);
sampleSet = sobolSet;
numSims = n;

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

%% Analysis
%start & end times for calculating steady state values
tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
tlong  = [285 300];     % for CO and SV determination

% Fitting PN_NA to VNS
VNStimepre = [120 180]; %after steady state reached, before VNS
VNStime     = [180 240]; %during VNS
VNStimepost = [240 300];

% Heart rates during different time periods
avgHRpre= zeros(1,n);
avgHRstim = zeros(1,n);
avgHRpost = zeros(1,n);

x = [1.59;1.54;1.62;1.59;1.54;1.62;1.59;1.54;1.62;67.5503;2.887;-0.0767;0.0018]; % tuning targets using poisson dist to estimate
% old tuning targets [1.59;1.57;1.62;1.59;1.57;1.62;1.59;1.57;1.62;67.5503;2.887];%84.6;2.392]; % Baseline, stimulation, recovery period averaged ICN neuron response to vagal stimulation, hR and Emaxlv from Ursino 1998, which park used for fitting
% Average across 7 neurons recorded in Rajendran data was used because all
% ICN neurons should roughly fire at this frequency
% Numbers are repeated three times, for NAPN,DMVPN, and LCN
% Park model output HR and Emaxlv: 67.5503;2.8871
simRes = zeros(length(x),n);
mse = zeros(1,n);

for i = 1:numSims
    simTime         = simout(1,i).time;
    tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
    tlongIdx 	= find(simTime >= tlong(1) & simTime <= tlong(2));
    tVNSpre 	= find(simTime >= VNStimepre(1) & simTime <= VNStimepre(2));
    tVNS        = find(simTime >= VNStime(1) & simTime <= VNStime(2));
    tVNSpost 	= find(simTime >= VNStimepost(1) & simTime <= VNStimepost(2));

    
    
    failedSim = isnan(avgHRpre);
    if failedSim == 1
        for j = 1:length(x)
            simRes(j,i) = NaN;
        end
    else
        simTime         = simout(1,i).time;
        tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
        tlongIdx        = find(simTime>= tlong(1)  & simTime <= tlong(2));
        
        Tperiod_vals    = simout(1,i).Tperiod;
        HR_vals         = simout(1,i).HR_TS;
        Emaxlv_vals     = simout(1,i).Emaxlv;
        
        Rsp_vals        = simout(1,i).Rsp;
        Rep_vals        = simout(1,i).Rep;
        Rmp_vals        = simout(1,i).Rmp;
        
        Vlv_vals        = simout(1,i).Vlv;
        
        Psa_vals        = simout(1,i).Psa;
        
        Flow_lv_vals    = simout(1,i).Flow_lv;
        Flow_rv_vals    = simout(1,i).Flow_rv;
        
        Pla_vals        = simout(1,i).Pla;
        Plv_vals        = simout(1,i).Plv;
        
        phi_vals        = simout(1,i).Phi;
        
        PN_NA_ff	= simout(1,i).PN_NA_ff;
        PN_DMV_ff	= simout(1,i).PN_DMV_ff;
        LCN_ff      = simout(1,i).LCN_ff;
        fevHR		= simout(1,i).fevHR;
        fevEmax		= simout(1,i).fevEmax;
        
        try
            CO_SV = fx_CO_SV(tlong, simTime, Flow_lv_vals, Flow_rv_vals, Tperiod_vals);
        catch
            %warning('error in CO/SV calculation likely due to parameters used');
            CO_SV = struct('RCO', NaN, 'RSV', NaN, 'LCO', NaN,  'LSV', NaN);
        end
        
        EDVthreshold    = 10; % end diastolic volume min threshold
        try
            [Vlv_pks, Vlv_index]    = findpeaks(Vlv_vals(tshortIdx), 'MinPeakHeight', EDVthreshold,  'MinPeakProminence', 5);   % returns peak values and indices of peaks
        catch
            warning('error in left ventricular volume calculation')
            Vlv_pks     = 999;
        end
        
        meanVlvPks      = mean(Vlv_pks);
        SysThreshold    = 10;
        
        try
            [Psys_pks, Psys_index]  = findpeaks(Psa_vals(tlongIdx),    'MinPeakHeight', SysThreshold, 'MinPeakProminence', 5);
            [Pdia_pks, Pdia_index]  = findpeaks(-1*Psa_vals(tlongIdx), 'MinPeakProminence', 5);
            
            if length(Psys_pks) ~= length(Pdia_pks)
                minLength   = min([length(Psys_pks), length(Pdia_pks)]);
                Psys_pks    = Psys_pks(1:minLength);
                Pdia_pks    = Pdia_pks(1:minLength);
            end
        catch
            Psys_pks    = 999*ones(10, 1);  % arbitrarily chosen vector length of 10
            Pdia_pks    = 999*ones(10, 1);  % arbitrarily chosen vector length of 10
        end
        
        % identify appropriate indices for Left Ventricualr End Systolic Volume
        % (ESV) and End Systolic Pressure (ESP)
        Vlv_vals_tshort     = Vlv_vals(tshortIdx);
        phi_vals_tshort     = phi_vals(tshortIdx);
        Vlv_idx             = find(Vlv_vals_tshort <= (meanVlvPks-10));
        phi_idx             = find(phi_vals_tshort <= 0.2);
        int_idx             = intersect(Vlv_idx, phi_idx);
        ESV_idx             = int_idx(find(diff(int_idx)> 5) + 1);
        ESV_vals            = Vlv_vals_tshort(ESV_idx);
        
        Plv_vals_tshort     = Plv_vals(tshortIdx);
        Vlv_idx2            = find(diff(Vlv_vals_tshort) <= 0.01 & diff(Vlv_vals_tshort) >= -0.01);
        int_idx2            = intersect(Vlv_idx,Vlv_idx2);
        ESP_idx             = int_idx2(find(diff(int_idx2) > 20)+1);
        ESP_vals            = Plv_vals_tshort(ESP_idx);
        
        % identify appropriate indices for Left Ventricualr End DIastolic Volume
        % (EDV) and End Diastolic Pressure (ESP)
        phi_idx             = find(phi_vals_tshort == 0);
        EDV_idx             = phi_idx(find(diff(phi_idx) > 2));
        EDV_vals            = Vlv_vals_tshort(EDV_idx);
        
        Vlv_idx3            = find(Vlv_vals_tshort > (meanVlvPks-10));
        int_idx3            = intersect(phi_idx, Vlv_idx3);
        EDP_idx             = int_idx3(find(diff(int_idx3)>2));
        EDP_vals            = Plv_vals_tshort(EDP_idx);
        
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
    disp('Parameters used: fmin, fmax, fmid, gain for NA, NActr, BR')
    sprintf('% 6.6f',sampleSet(minidx,:)')
    disp('Initial HR')
%     avgHRbase(minidx)
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

%% display HR, CO, ejection fraction, etc
Params = paramMatrix(1,:);
simIn = Simulink.SimulationInput(mdlName);
    simIn = simIn.setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);
    
    simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...        % Set ICN parameters
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
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(PNDMVdelay), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(Params(15)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)));
    
    simout = sim(simIn);

disp('Systolic P, Diastolic P, HR, CO, EjFr')
PhysOutputs_Gen(simout)


%% Exit code
close_system
delete(myPool);
exit