function [mat_filename, HRvalidation, MAP_EI_mdl, SP_store, DP_store, CO_store, ...
    NA_activity_postIR, DMV_activity_postIR, symp_activity_postIR, ...
    SDRR_store, RMSSD_store, pNN50_store] = baroreflex_curve(mdlName, ...
    BPvalues,filename,ICNparams,NAparams,NTSparams,heart_params,ka,baseline)
% Baroreflex curve plotted as a function of changing Psa and RR interval
% response
% Inputs:
% mdlName: string with name of model to be simulated
% BPvalues: row vector of neck chamber pressures. carotid distending
% pressure = systolic pressure - neck chamber pressure (Rea1987)
% plotName: string with name for baroreflex curve plot leaving out the file suffix. ei 'plot'
% ICNparams: 1x17 vector of parameters
% NAparams: 1x12 vector of parameters
% NTSparams: 1x12 vector of parameters
% heart_params: 1x5 vector
% ka: 1x1 vector
% baseline: if baseline = 1, calculate baseline metrics. if baseline = 0,
% don't calculate baseline metrics

%% Simulation input
BPswitch = 1; % use BP set
kRSA = 0.5;

numParamSets = length(ka);
numBPvals = length(BPvalues);
numSims = numBPvals*numParamSets;

simIn(1:numBPvals) = Simulink.SimulationInput(mdlName);
simIn = repmat(simIn,numParamSets,1);
for j = 1:numParamSets

    for idx = 1:numBPvals
        simIn(j,idx) = simIn(j,idx).setModelParameter('SaveTime', 'on', ...
            'SaveOutput', 'on', ...
            'TimeOut', 240);
        simIn(j,idx) = simIn(j,idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(ICNparams(j,1)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(ICNparams(j,2)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(ICNparams(j,3)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(ICNparams(j,4)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(ICNparams(j,5)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(ICNparams(j,6)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(ICNparams(j,7)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(ICNparams(j,8)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(ICNparams(j,9)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(ICNparams(j,10)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(ICNparams(j,11)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(ICNparams(j,12)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(ICNparams(j,13)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(ICNparams(j,14)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(ICNparams(j,15)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(ICNparams(j,16)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(ICNparams(j,17)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(kRSA), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmin'], 'Value', num2str(NTSparams(j,1)), ...         % NTSparams
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmax'], 'Value', num2str(NTSparams(j,2)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmid'], 'Value', num2str(NTSparams(j,3)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/k, gain'], 'Value', num2str(NTSparams(j,4)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmin'], 'Value', num2str(NTSparams(j,5)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmax'], 'Value', num2str(NTSparams(j,6)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmid'], 'Value', num2str(NTSparams(j,7)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/k, gain'], 'Value', num2str(NTSparams(j,8)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmin'], 'Value', num2str(NTSparams(j,9)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmax'], 'Value', num2str(NTSparams(j,10)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmid'], 'Value', num2str(NTSparams(j,11)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/k, gain'], 'Value', num2str(NTSparams(j,12)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmin'], 'Value', num2str(NAparams(j,1)), ...     % NAparams
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmax'], 'Value', num2str(NAparams(j,2)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmid'], 'Value', num2str(NAparams(j,3)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/k, gain'], 'Value', num2str(NAparams(j,4)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmin'], 'Value', num2str(NAparams(j,5)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmax'], 'Value', num2str(NAparams(j,6)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmid'], 'Value', num2str(NAparams(j,7)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/k, gain'], 'Value', num2str(NAparams(j,8)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmin'], 'Value', num2str(NAparams(j,9)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmax'], 'Value', num2str(NAparams(j,10)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmid'], 'Value', num2str(NAparams(j,11)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/k, gain'], 'Value', num2str(NAparams(j,12)), ...
            [mdlName '/Left Heart/Cla'], 'Value', num2str(heart_params(j,1)), ...                                    % heart_params
            [mdlName '/Autonomic Nervous System/Left ventricle/P0lv'], 'Value', num2str(heart_params(j,2)), ...
            [mdlName '/Autonomic Nervous System/Left ventricle/kElv'], 'Value', num2str(heart_params(j,3)), ...
            [mdlName '/Autonomic Nervous System/Left ventricle/Vulv'], 'Value', num2str(heart_params(j,4)), ...
            [mdlName '/Left Heart/Subsystem2/kRlv'], 'Value', num2str(heart_params(j,5)), ...
            [mdlName '/Autonomic Nervous System/Emaxlv_inv/Emaxlv0_denervated (Emaxlv0)'], 'Value', num2str(heart_params(j,6)), ...
            [mdlName '/Autonomic Nervous System/Symp. Efferent Pathways/fes_inf'], 'Value', num2str(heart_params(j,7)), ...
            [mdlName '/Autonomic Nervous System/Carotid Sinus (baroreceptors)/ka'], 'Value', num2str(ka(j)), ...     % ka
            [mdlName '/Autonomic Nervous System/BP'], 'Value', num2str(BPvalues(idx)), ...
            [mdlName '/Autonomic Nervous System/BPswitch'], 'Value', num2str(BPswitch));

    end
end

%% Run simulations
% simOut = parsim(simIn);
simOut = sim(simIn);

%% Analysis
% determine % HR and BP change for each model
HR = zeros(numParamSets,numBPvals);

tshort = [150 155]; % BP step at 150 s
for i = 1:numBPvals
    for j = 1:numParamSets
        try
            HR(j,i) = mean(getsampleusingtime(simOut(j,i).HR_TS,tshort(1),tshort(2)),'Weighting','time');
        catch
            disp('Warning: simulation failed. HR value assigned as nan')
            HR(j,i) = nan;
        end
    end
end

% Calculate delta RR interval
% Baseline is 120 mm Hg RR interval
idx_base = find(BPvalues == 0);
RR = 60./HR; % convert to RR interval in s

% Calculate RR interval change from RR interval at 0 neck chamber pressure
RR_120 = RR(idx_base); % RR interval at 0 neck chamber pressure
deltaRR = 1000*(RR - RR_120); % change in RR interval from 0 neck chamber pressure in ms
if length(deltaRR) == 0
    disp('Error: BPvalues must include 0 mm Hg neck chamber pressure as baseline for comparison')
end

% Calculate systolic blood
tshort = [150 155];
tlong = [150 155];
SBP = zeros(numParamSets,numBPvals);
DBP = zeros(numParamSets,numBPvals);
MAP = zeros(numParamSets,numBPvals);
for i = 1:numBPvals
    for j = 1:numParamSets
        simRes = PhysOutputs_Gen_TS(simOut(j,i),tshort,tlong);
        SBP(j,i) = simRes(1);
        DBP(j,i) = simRes(2);
        MAP(j,i) = DBP(1,i) + (1/3)*(SBP(1,i)-DBP(1,i));
    end
end
CSP = SBP + BPvalues; % carotid sinus pressure
baselineMAP = MAP(idx_base);


%% Optional calculation of baseline hemodynamic metrics
if baseline == 1
    %start & end times for calculating steady state values
    tshort = [150 165]; % chose 15 sec range since this is what is used clinically to determine BPM
    tlong = [150 180]; % for CO and SV determination

    % Preallocate
    HRvalidation = zeros(1,numParamSets);
    MAP_EI_mdl = zeros(1,numParamSets);
    NA_activity_postIR = zeros(1,numParamSets);
    DMV_activity_postIR = zeros(1,numParamSets);
    symp_activity_postIR = zeros(1,numParamSets);
    SP_store = zeros(1,numParamSets);
    DP_store = zeros(1,numParamSets);
    CO_store = zeros(1,numParamSets);
    SDRR_store = zeros(1,numParamSets);
    RMSSD_store = zeros(1,numParamSets);
    pNN50_store = zeros(1,numParamSets);

    zeroIdx = find(BPvalues == 0);
    % Calculate HR and MAP for each parameter set
    for i = 1:numParamSets
        if ~isempty(simOut(i,zeroIdx).ErrorMessage)
            HRvalidation(i) = nan;
            MAP_EI_mdl(i) = nan;
            SP_store(i) = nan;
            DP_store(i) = nan;
            NA_activity_postIR(i) = nan;
            DMV_activity_postIR(i) = nan;
            % CHANGED
            symp_activity_postIR(i) = nan;

        else
            HR_vals         = simOut(i,zeroIdx).HR_TS;
            avgHRpre = mean(getsampleusingtime(HR_vals,165,180),'Weighting','time');


            failedSim = isnan(avgHRpre);
            if failedSim == 1
                HRvalidation(i) = nan;
            else

                % Calculate percent change in heart rate due to step
                if length(mean(getsampleusingtime(HR_vals,165,180),'Weighting','time')) == 0
                    HRvalidation(i)      = nan; %baseline
                    MAP_EI_mdl(i)      = nan; %baseline
                else
                    physOutputs = PhysOutputs_Gen_TS(simOut(i,zeroIdx),tshort,tlong);
                    SP = physOutputs(1);
                    DP = physOutputs(2);
                    CO = physOutputs(4);
                    HRvalidation(i) = mean(getsampleusingtime(HR_vals,150,180),'Weighting','time');
                    MAP_EI_mdl(i) = (SP + 2*DP)/3;
                    SP_store(i) = SP;
                    DP_store(i) = DP;
                    CO_store(i) = CO;
                    NA_activity_postIR(i) = mean(getsampleusingtime(simOut(i,zeroIdx).fevHR,tshort(1),tshort(2)),'Weighting','time');
                    DMV_activity_postIR(i) = mean(getsampleusingtime(simOut(i,zeroIdx).fevEmax,tshort(1),tshort(2)),'Weighting','time');
                    symp_activity_postIR(i) = mean(getsampleusingtime(simOut(i,zeroIdx).fesh,tshort(1),tshort(2)),'Weighting','time');

                    trange = [120 200];
                    [SDRR, RMSSD, pNN50,~,~,~,~,~,~] = HRVmetrics(simOut(i,zeroIdx),trange);
                    SDRR_store(1,i) = SDRR;
                    RMSSD_store(1,i) = RMSSD;
                    pNN50_store(1,i) = pNN50;
                end
            end
        end
    end
else
    HRvalidation= [];
    MAP_EI_mdl= [];
    SP_store = [];
    DP_store = [];
    CO_store = [];
    NA_activity_postIR = [];
    DMV_activity_postIR = [];
    symp_activity_postIR = [];

    SDRR_store = [];
    RMSSD_store = [];
    pNN50_store = [];
end
%% plot baroreflex curve and save data
% plot formatting
% fs = 32;
% ms = 12;

% Osculati 1990 mean data
% load data
% neckChamberPressure = [8	13.2	18.5	23.7	30	37.5]; % mm hg
% RRbase = 777.8; % +/- 34.7 msec
% deltaRR = [46.50360279	64.79801636	89.04640376	126.5452133	144.0065953	171.9898378];
% SEM = [7.550265772	10.06205045	13.11093148	24.29514847	24.91603345	24.7516633]; % msec
% num = 12;
% SD = SEM.* sqrt(num); % standard deviation
% CI95 = 1.96.* SD;
% RR_Osculati = RRbase + deltaRR;

%% Plot
% data
% x = neckChamberPressure;
% y = RR_Osculati;
% errors = CI95;

% figure;
% hold on;
% plot(x, y, 'r', 'LineWidth', 2);
%
% % Calculate upper and lower bounds for the shaded region
% upper_bound = y + errors;
% lower_bound = y - errors;
%
% % Plot the shaded region
% x_fill = [x, fliplr(x)];
% y_fill = [upper_bound, fliplr(lower_bound)];
% fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
%
% % Add labels and title
% xlabel('Neck chamber pressure (mm Hg)');
% ylabel('RR interval (ms)');
% % title('Line with Shaded Region (Standard Error)');
% purple = [0.4940 0.1840 0.5560];
% plot(BPvalues,1000.*60./HR,'b-','MarkerSize',ms,'MarkerFaceColor','b')
% set(gca,'FontSize',fs)
%
% set(gcf, 'Position',  [10, 10, 600, 500])
% file_suffix = '.png';
% png_filename = [filename file_suffix];
% saveas(gcf,png_filename)

% save data
file_suffix = '.mat';
mat_filename = [filename file_suffix];
save(mat_filename,'SBP','deltaRR','BPvalues','HR','baselineMAP') % HR a column vector
end
