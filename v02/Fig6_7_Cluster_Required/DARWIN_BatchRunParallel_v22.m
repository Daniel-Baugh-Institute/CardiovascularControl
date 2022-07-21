% Timeseries simOut format for averaging ICN ff, HR, and Emaxlv

clear;
%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);


%% Define name of simulink model
mdlName     = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v11';

%% Create parameter set and number of simulations
n = 5000; % number of sobol sample sets
[sobolSet] = SobolSample([0.06,1.3,2,0.1,0.06,1.3,2,0.1,0.06,5,5,0.01,1,1,0.001], [10,150,22,10,25,150,300,5,5,21,15,4,20,20,0.999],15,n);
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
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,15)));
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

x = [1.59;1.54;1.62;1.59;1.54;1.62;1.59;1.54;1.62;67.5503;2.887];% tuning targets using poisson dist to estimate
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

%     avgHRbase(i) = mean(getsampleusingtime(simout(1,i).HR,165,180),'Weighting','time');
%     avgHRstim(i) = mean(getsampleusingtime(simout(1,i).HR,240),'Weighting','time');
%     avgHRrecov(i) = mean(getsampleusingtime(simout(1,i).HR,285,300),'Weighting','time');
    
    
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
        HR_vals         = simout(1,i).HR;
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
        if length(mean(getsampleusingtime(PN_NA_ff,240,300),'Weighting','time')) == 0
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
[sobolSet] = SobolSample(minParam', maxParam', 15, n);
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
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(sampleSet(idx,15)));
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

x = [1.59;1.54;1.62;1.59;1.54;1.62;1.59;1.54;1.62;67.5503;2.887]; % tuning targets using poisson dist to estimate
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

%     avgHRbase(i) = mean(getsampleusingtime(simout(1,i).HR,165,180),'Weighting','time');
%     avgHRstim(i) = mean(getsampleusingtime(simout(1,i).HR,240),'Weighting','time');
%     avgHRrecov(i) = mean(getsampleusingtime(simout(1,i).HR,285,300),'Weighting','time');
    
    
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
        HR_vals         = simout(1,i).HR;
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
        if length(mean(getsampleusingtime(PN_NA_ff,240,300),'Weighting','time')) == 0
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
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(15)));
    
    simout = sim(simIn);

disp('Systolic P, Diastolic P, HR, CO, EjFr')
PhysOutputs_Gen(simout)

%% subpanel 2: RR interval and emaxlv
simOut = simout;
tplot = [165 180];
simTime         = simOut.time;
tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));

% RR interval calculation
time = simOut.time(tplotIdx);
phi = simOut.Phi(tplotIdx);

% Determine indices for when phi=0 (beginning and end of heart beat)
for i = length(phi):-1:2
    if phi(i) == phi(i-1)
        phi(i-1) = [];
        time(i-1) = [];
    end
end

RRidx = find(~phi);
RRtimes = time(RRidx);
RRinterval = diff(RRtimes);


% Plot as staircase
timePlot = repelem(RRtimes,2);
timePlot(1) = []; % get rid of first entry and last time entry to align 
timePlot(end) = [];
RRintervalPlot = repelem(RRinterval,2);

% Plot
plotlim = [165 180];
fs = 18;
lw = 2;

fig = figure(1);
subplot(2,1,1)
plot(timePlot,RRintervalPlot,'LineWidth',2)
xlabel('Time (s)','FontSize',fs)
ylabel('RR Interval (s)','FontSize',fs)
ax = gca;
set(gca,'FontSize',fs)
%title('Heart Rate','FontSize',18)
xlim(plotlim)

subplot(2,1,2)
plot(getsampleusingtime(Emaxlv_vals,165,180),'LineWidth',lw)
xlabel('Time (s)','FontSize',fs)
ylabel('Emax_{LV} (mmHg/mL)','FontSize',fs)
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

set(gcf, 'Position',  [100, 100, 600, 500])

% Make folder with timestamp
format longG
t = now
timestamp = num2str(t);
command = ['mkdir ' timestamp];
status = system(command);
cd(timestamp) 
d = pwd;
% Save file
saveas(gcf,fullfile(d,'plot_ICNhr'),'png');


%% subpanel 1: ICN input-output firing frequencies (Base and no local circuit)
% NA-PN firing frequency
figure(2);

% Output firing frequencies
subplot(3,2,2)
plot(getsampleusingtime(PN_NA_ff,165,180),'LineWidth',lw)
title({'Output Firing Frequencies'},'FontSize',fs)
set(gca,'FontSize',fs)
xlim(plotlim)

% DMV-PN firing frequency
subplot(3,2,4)
plot(getsampleusingtime(PN_DMV_ff,165,180),'LineWidth',lw)
set(gca,'FontSize',fs)
xlim(plotlim)

%LCN firing frequency
subplot(3,2,6)
plot(getsampleusingtime(LCN_ff,165,180),'LineWidth',lw)
xlabel('Time (s)','FontSize',fs)
set(gca,'FontSize',fs)
xlim(plotlim)

% Input firing frequencies
subplot(3,2,1)
plot(simOut.time(tplotIdx),simOut.fevHR(tplotIdx),'LineWidth',lw)
ylabel({'Firing'; 'Frequency (Hz)'},'FontSize',fs)
title({'Input Firing Frequencies'},'FontSize',fs)
xlim(plotlim)
set(gca,'FontSize',fs)

subplot(3,2,3)
plot(simOut.time(tplotIdx),simOut.fevEmax(tplotIdx),'LineWidth',lw)
ylabel({'Firing'; 'Frequency (Hz)'},'FontSize',fs)
set(gca,'FontSize',fs)
xlim(plotlim)

subplot(3,2,5)
plot(simOut.time(tplotIdx),simOut.LCNin(tplotIdx),'LineWidth',lw)
xlabel('Time (s)','FontSize',fs)
ylabel({'Net'; 'Activity (Hz)'},'FontSize',fs)
set(gca,'FontSize',fs)
xlim(plotlim)

set(gcf, 'Position',  [100, 100, 1200, 600])
saveas(gcf,fullfile(d,'plot_ICNff'),'png');

%% PV plot
figure(9)
plotlim = [0 200];
plot(simOut.Vlv(tplotIdx),simOut.Plv(tplotIdx),'b-','LineWidth',lw)
xlabel('Left Ventricular Volume (mL)','FontSize',fs)
ylabel('Left Ventricular Pressure (mmHg)','FontSize',fs)
set(gcf, 'Position',  [0, 0, 600, 600])
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)
saveas(gcf,fullfile(d,'plot_PV_overlay'),'png');
ylim(plotlim)

% %% Plots for Suga 1976 and Greenwood 1980 Data
% % Return to previous folder with model files
% cd ..
% 
% mdlName = 'ICN_with_BR_input_model4test_clusterParams_PF2C_2lane_v2';
% % HR
% % Vary input arterial pressure in mmHg
% Psa = 50:25:150;
% 
% HR2lanes = zeros(length(Psa));
% Emax2lanes = zeros(length(Psa));
% 
% 
% % Simulation input with tuned model parameters and input Psa
% simIn = Simulink.SimulationInput(mdlName);
% 
% % Two lanes
% simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(Params(2)), ... %output is the fmax value %(deltaNA), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(Params(3)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(Params(4)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(Params(5)), ... % NActr
%     [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(Params(6)), ... %deltaNActr), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(Params(7)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(Params(8)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(Params(9)), ... %DMV
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(Params(10)), ... %(deltaDMV), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(Params(11)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(Params(12)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(Params(13)), ... % BR cell grp
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(14)), ... %(deltaBR), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(15)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)));
% 
% % Run simulation
% simOutBP = sim(simIn);
% 
% % Indices for each arterial pressure
% t50 = [105 120];     % chose 15 sec range since this is what is used clinically to determine BPM
% t75 = [225 240];
% t100 = [345 360];
% t125 = [465 480];
% t150 = [585 600];
% 
% 
% 
% simTime = simOutBP.time;
% 
% t50Idx          = find(simTime>= t50(1) & simTime <= t50(2));
% t75Idx          = find(simTime>= t75(1) & simTime <= t75(2));
% t100Idx         = find(simTime>= t100(1) & simTime <= t100(2));
% t125Idx          = find(simTime>= t125(1) & simTime <= t125(2));
% t150Idx          = find(simTime>= t150(1) & simTime <= t150(2));
% 
% 
% HR2lanes(1) = mean(simOutBP.HR(t50Idx));
% HR2lanes(2) = mean(simOutBP.HR(t75Idx));
% HR2lanes(3) = mean(simOutBP.HR(t100Idx));
% HR2lanes(4) = mean(simOutBP.HR(t125Idx));
% HR2lanes(5) = mean(simOutBP.HR(t150Idx));
% 
% 
% Emax2lanes(1) = mean(simOutBP.Emaxlv(t50Idx));
% Emax2lanes(2) = mean(simOutBP.Emaxlv(t75Idx));
% Emax2lanes(3) = mean(simOutBP.Emaxlv(t100Idx));
% Emax2lanes(4) = mean(simOutBP.Emaxlv(t125Idx));
% Emax2lanes(5) = mean(simOutBP.Emaxlv(t150Idx));
% 
% 
% %% % No Local reflex
% ParamsNLR = Params;
% ParamsNLR(15) = [];
% 
% mdlName = 'ICN_with_BR_input_model4test_clusterParams_PF2C_NLR';
% 
% % Simulation input with tuned model parameters and input Psa
% simIn = Simulink.SimulationInput(mdlName);
% 
% simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(ParamsNLR(1)), ...        % Set ICN parameters
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(ParamsNLR(2)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(ParamsNLR(3)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(ParamsNLR(4)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(ParamsNLR(5)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(ParamsNLR(6)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(ParamsNLR(7)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(ParamsNLR(8)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(ParamsNLR(9)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(ParamsNLR(10)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(ParamsNLR(11)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(ParamsNLR(12)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(ParamsNLR(13)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(ParamsNLR(14)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(ParamsNLR(15)));
% 
% % Run simulation
% simOut = sim(simIn);
% 
% 
% % Indices for each arterial pressure
% simTime = simOut.time;
% 
% t50Idx          = find(simTime>= t50(1) & simTime <= t50(2));
% t75Idx          = find(simTime>= t75(1) & simTime <= t75(2));
% t100Idx         = find(simTime>= t100(1) & simTime <= t100(2));
% t125Idx          = find(simTime>= t125(1) & simTime <= t125(2));
% t150Idx          = find(simTime>= t150(1) & simTime <= t150(2));
% 
% 
% HRnoLCN(1) = mean(simOut.HR(t50Idx));
% HRnoLCN(2) = mean(simOut.HR(t75Idx));
% HRnoLCN(3) = mean(simOut.HR(t100Idx));
% HRnoLCN(4) = mean(simOut.HR(t125Idx));
% HRnoLCN(5) = mean(simOut.HR(t150Idx));
% EmaxnoLCN(1) = mean(simOut.Emaxlv(t50Idx));
% EmaxnoLCN(2) = mean(simOut.Emaxlv(t75Idx));
% EmaxnoLCN(3) = mean(simOut.Emaxlv(t100Idx));
% EmaxnoLCN(4) = mean(simOut.Emaxlv(t125Idx));
% EmaxnoLCN(5) = mean(simOut.Emaxlv(t150Idx));
% 
% %% Data from Suga 1976
% HRpercent = [104.73492723492723 106.23180873180868 99.39362439362435 90.63582813582812 75.46257796257794]/ 100;
% Emaxpercent = [112.18222680716232 109.95320904870681 99.39577039274926 85.84113182521554 76.9471667526343] / 100;
% BaseHR = 61.2; % Used by Park based on Ursino
% BaseEmax = 2.695; % Used by Park based on Ursino
% HRSuga = BaseHR * HRpercent;
% EmaxSuga = BaseEmax * Emaxpercent;
% 
% 
% %% Plots
% fig = figure(1);
% ms = 6;
% lw = 3;
% fs = 16;
% subplot(1,2,1)
% plot(Psa,HRSuga,'ko','MarkerSize',ms,'LineWidth',lw)
% hold on
% plot(Psa,HR2lanes(:,1),'b-','MarkerSize',ms,'LineWidth',lw)
% hold on
% plot(Psa,HRnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
% ylabel('HR (bpm)','FontSize',fs)
% ylim([40 120])
% ax = gca;
% set(gca,'FontSize',fs)
% 
% 
% 
% subplot(1,2,2)
% plot(Psa,EmaxSuga,'ko','MarkerSize',ms,'LineWidth',lw)
% hold on
% plot(Psa,Emax2lanes(:,1),'b-','MarkerSize',8,'LineWidth',lw)
% hold on
% plot(Psa,EmaxnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
% ylabel('E_{max,lv} (mmHg/mL)','FontSize',fs)
% %legend('Corresponding Experimental Data','Two Lanes','No Local Reflex','Location','northeastoutside')
% %ylim([0.5 4.5])
% ax = gca;
% set(gca,'FontSize',fs)
% %ax.FontSize = 14;
% 
% % common x axis label
% han = axes(fig,'visible','off');
% han.XLabel.Visible = 'on';
% %xlabel(han,'Mean Arterial Pressure (mmHg)','FontSize',fs)
% set(gcf, 'Position',  [50, 50, 1000, 400])
% saveas(gcf,fullfile(d,'ParkF2C_AP'),'png');
% hold off
% %% Vary Lung Tidal Volume
% mdlName = 'ICN_with_BR_input_model4test_clusterParams_ParkFig2C_Vlung';
% Vlung = [0, 0.75, 1.25, 1.75]; % Percent increase from nominal
% 
% % Park model params
% BRclosedLoop_model_initialize_params;
% 
% 
% HR2lanes = zeros(length(Vlung));
% Emax2lanes = zeros(length(Vlung));
% 
% 
% % Simulation input with tuned model parameters and input Psa
% simIn = Simulink.SimulationInput(mdlName);
% 
% simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(Params(2)), ... %output is the fmax value %(deltaNA), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(Params(3)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(Params(4)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(Params(5)), ... % NActr
%     [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(Params(6)), ... %deltaNActr), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(Params(7)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(Params(8)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(Params(9)), ... %DMV
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(Params(10)), ... %(deltaDMV), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(Params(11)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(Params(12)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(Params(13)), ... % BR cell grp
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(14)), ... %(deltaBR), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(15)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)));
% 
% % Run simulation
% simOut = sim(simIn);
% 
% % Indices for each arterial pressure
% t50 = [105 120];     % chose 15 sec range since this is what is used clinically to determine BPM
% t75 = [225 240];
% t100 = [345 360];
% t125 = [465 480];
% 
% 
% simTime = simOut.time;
% 
% t50Idx          = find(simTime>= t50(1) & simTime <= t50(2));
% t75Idx          = find(simTime>= t75(1) & simTime <= t75(2));
% t100Idx         = find(simTime>= t100(1) & simTime <= t100(2));
% t125Idx          = find(simTime>= t125(1) & simTime <= t125(2));
% 
% 
% 
% HR2lanes(1) = mean(simOut.HR(t50Idx));
% HR2lanes(2) = mean(simOut.HR(t75Idx));
% HR2lanes(3) = mean(simOut.HR(t100Idx));
% HR2lanes(4) = mean(simOut.HR(t125Idx));
% 
% Emax2lanes(1) = mean(simOut.Emaxlv(t50Idx));
% Emax2lanes(2) = mean(simOut.Emaxlv(t75Idx));
% Emax2lanes(3) = mean(simOut.Emaxlv(t100Idx));
% Emax2lanes(4) = mean(simOut.Emaxlv(t125Idx));
% 
% %% % No Local reflex
% mdlName = 'ICN_with_BR_input_model4test_clusterParams_ParkFig2C_VlungNLR';%'ICN_with_BR_input_model4test_clusterParams_ICNnoLCN_ParkF2_lung';
% 
% % Simulation input with tuned model parameters and input Psa
% simIn = Simulink.SimulationInput(mdlName);
% 
% simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(ParamsNLR(1)), ...        % Set ICN parameters
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(ParamsNLR(2)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(ParamsNLR(3)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(ParamsNLR(4)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(ParamsNLR(5)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(ParamsNLR(6)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(ParamsNLR(7)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(ParamsNLR(8)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(ParamsNLR(9)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(ParamsNLR(10)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(ParamsNLR(11)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(ParamsNLR(12)), ...
%     [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(ParamsNLR(13)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(ParamsNLR(14)), ...
%     [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(ParamsNLR(15)));
% 
% % Run simulation
% simOut = sim(simIn);
% 
% % Indices for each arterial pressure
% 
% simTime = simOut.time;
% 
% t50Idx          = find(simTime>= t50(1) & simTime <= t50(2));
% t75Idx          = find(simTime>= t75(1) & simTime <= t75(2));
% t100Idx         = find(simTime>= t100(1) & simTime <= t100(2));
% t125Idx          = find(simTime>= t125(1) & simTime <= t125(2));
% 
% 
% 
% HRnoLCN(1) = mean(simOut.HR(t50Idx));
% HRnoLCN(2) = mean(simOut.HR(t75Idx));
% HRnoLCN(3) = mean(simOut.HR(t100Idx));
% HRnoLCN(4) = mean(simOut.HR(t125Idx));
% 
% EmaxnoLCN(1) = mean(simOut.Emaxlv(t50Idx));
% EmaxnoLCN(2) = mean(simOut.Emaxlv(t75Idx));
% EmaxnoLCN(3) = mean(simOut.Emaxlv(t100Idx));
% EmaxnoLCN(4) = mean(simOut.Emaxlv(t125Idx));
% 
% %%
% % Data from Suga 1976
% % HRpercent = [104.73492723492723 106.23180873180868 99.39362439362435 90.63582813582812 75.46257796257794]/ 100;
% % Emaxpercent = [112.18222680716232 109.95320904870681 99.39577039274926 85.84113182521554 76.9471667526343] / 100;
% % BaseHR = 61.2; % Used by Park based on Ursino
% % BaseEmax = 2.695; % Used by Park based on Ursino
% HRGreenwood = [68.67647058823529 77.6470588235294  82.94117647058822 90.14705882352939];
% HRPark = [76.02847548690396 80.21383478844861 84.92975151108124 93.75795836131628];
% xGreenwood = [0, 0.73015873015872930, 1.25396825396825480, 1.78571428571428650];
% xEmax = [0, 0.75, 1.25, 1.75];
% EmaxPark = [2.3947848045674136 2.2945652173913045 2.3663043478260875 2.378072024593764];
% EmaxGreenwood = [ 2.2539174352217834 2.0771739130434783 2.1228238910847606 2.173722002635046];
% %%
% % Plots
% ms = 6;
% lw = 3;
% fs = 16;
% 
% fig = figure(2);
% subplot(1,2,1)
% plot(xGreenwood,HRGreenwood,'ko','MarkerSize',ms,'LineWidth',lw)
% hold on
% % plot(xGreenwood,HRnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
% hold on
% plot(Vlung,HR2lanes(:,1),'b-','MarkerSize',ms,'LineWidth',lw)
% ylabel('HR (bpm)','FontSize',fs)
% ylim([40 120])
% ax = gca;
% set(gca,'FontSize',fs)
% %ax.FontSize = 14;
% 
% subplot(1,2,2)
% plot(xEmax,EmaxGreenwood,'ko','MarkerSize',ms,'LineWidth',lw)
% hold on
% % plot(xEmax,EmaxnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
% hold on
% plot(Vlung,Emax2lanes(:,1),'b-','MarkerSize',ms,'LineWidth',lw)
% ylabel('E_{max,lv} (mmHg/mL)','FontSize',fs)
% %ylim([0.5 4.5])
% ax = gca;
% set(gca,'FontSize',fs)
% %ax.FontSize = 14;
% 
% % common x axis label
% han = axes(fig,'visible','off');
% han.XLabel.Visible = 'on';
% %xlabel(han,'Lung Tidal Volume (% from nominal)','FontSize',fs)
% set(gcf, 'Position',  [50, 50, 1000, 400])
% hold off
% saveas(gcf,fullfile(d,'ParkFig2C_Vlung'),'png');
%% Exit code
close_system
delete(myPool);
exit