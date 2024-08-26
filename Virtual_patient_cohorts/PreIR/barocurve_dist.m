function [sampleSet, SP, DP, HR, MAP,A1Store,A2Store,A3Store,A4Store,NAactivity,DMVactivity,sympActivity] = barocurve_dist(BPvalues,parameters,filename)
% Given a parameter value or set of parameter values, what is the
% distribution of baroreflex curves? How does it compare to the data from
% individuals in Seredynski?

% Seredynski
% doi: 10.3389/fphys.2021.703692

% Michelle Gee
% February 23, 2024

% Function with inputs of patient
% Sample to select parameter sets
% For each parameter set, define simIn for a set of neck chamber pressures
% (do this instead of calling the baroreflex fuction so that it can be
% parallelized)
% Calculate the heart rate and MAP for the 5 second neck chamber
% application
% Store in matrix and save as .mat file

%%
% addpath('C:\Users\mmgee\MATLAB\Projects\untitled')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/Fig_2_Rea_baroreflex/')
addpath('C:\Users\mmgee\MATLAB\Projects\untitled')
my_dir = pwd;
addpath(genpath(my_dir))


%% Define name of simulink model
mdlName     = 'ICN_model_v15_Mastitskaya2012_control_r2020b';

%% Create parameter set and number of simulations
n = 500; % number of sobol sample sets


% Parameter values that significantly affect HR and/or MAP from sensitivity analysis by neural group:
% CP_fmax (23), BR_fmid (20), CP_fmid (21), LS_fmid (28), LS_fmax (27),
% PNDMV_fmin (9), PNNA_fmid (3), PNNA_fmin (1), PNNA_fmax (2), PNDMV_fmid (11), PNNA_k (4)
% NA_fmin (30), NA_fmid (32), DMV_fmid (40), ka (47), DMV_fmax (39), NActr_fmid
% ka
% Cla, P0lv, kElv, Vulv, fesinf

% Indices of parameter values from sensitivity analysis by neural group:
sigParamIdx = [1, 2, 3, 4, 9, 11, 20, 21, 23, 27, 28, 30, 32, 39, 40, 42, 43, 44, 45, 46, 47];

% Significant parameters based on combined sensitivity analysis:
sigParamCombinedIdx = [1 2 3 4 5 7 8 9 11 16 20 23 27 28 30 31 32 36 39 40 42 45 46 47 48 49 50 51];
%42-48 are
% P0lv = 1.5;
% kElv = 0.014;
% Vulv = 16.77;
% kRlv = 3.75e-4
% Emaxlv0 = 1.283;
% fes_inf = 2.1;

kRSA = 0.5;
paramVec = parameters(sigParamCombinedIdx);
BPswitch = 1; % use open loop


% Set parameter range for Sobol sampling
factor = 2.5;
minParams = paramVec./factor;
maxParams = paramVec.*factor;
[sobolSet] = SobolSample(minParams,maxParams,length(paramVec),n);

sampleSet = sobolSet;
n = 500;
m = length(BPvalues);
zeroIdx = find(BPvalues == 0);
%% Create an array of simulation input objects and specify the sweep value for each simulation
simIn(1:m) = Simulink.SimulationInput(mdlName);
simIn = repmat(simIn, n,1);

for i = 1:n
    for j = 1:m
        simIn(i,j) = simIn(i,j).setModelParameter('SaveTime', 'on', ...
            'SaveOutput', 'on', ...
            'TimeOut', 240);

        simIn(i,j) = simIn(i,j).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(sampleSet(i,1)), ... % ICN params
            [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(sampleSet(i,2)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(sampleSet(i,3)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(sampleSet(i,4)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(sampleSet(i,5)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(parameters(6)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(sampleSet(i,6)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(sampleSet(i,7)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(sampleSet(i,8)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(parameters(10)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(sampleSet(i,9)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(parameters(12)), ...
            [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(parameters(13)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(parameters(14)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(parameters(15)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(sampleSet(i,10)), ...
            [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(parameters(17)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(kRSA), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmin'], 'Value', num2str(parameters(18)), ...         % NTSparams
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmax'], 'Value', num2str(parameters(19)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmid'], 'Value', num2str(sampleSet(i,11)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/k, gain'], 'Value', num2str(parameters(21)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmin'], 'Value', num2str(parameters(22)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmax'], 'Value', num2str(sampleSet(i,12)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmid'], 'Value', num2str(parameters(24)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/k, gain'], 'Value', num2str(parameters(25)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmin'], 'Value', num2str(parameters(26)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmax'], 'Value', num2str(sampleSet(i,13)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmid'], 'Value', num2str(sampleSet(i,14)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/k, gain'], 'Value', num2str(parameters(29)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmin'], 'Value', num2str(sampleSet(i,15)), ...     % NA/DMV params
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmax'], 'Value', num2str(sampleSet(i,16)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmid'], 'Value', num2str(sampleSet(i,17)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/k, gain'], 'Value', num2str(parameters(33)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmin'], 'Value', num2str(parameters(34)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmax'], 'Value', num2str(parameters(35)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmid'], 'Value', num2str(sampleSet(i,18)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/k, gain'], 'Value', num2str(parameters(37)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmin'], 'Value', num2str(parameters(38)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmax'], 'Value', num2str(sampleSet(i,19)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmid'], 'Value', num2str(sampleSet(i,20)), ...
            [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/k, gain'], 'Value', num2str(parameters(41)), ...
            [mdlName '/Left Heart/Cla'], 'Value', num2str(sampleSet(i,21)), ...
            [mdlName '/Autonomic Nervous System/Vula'], 'Value', num2str(parameters(43)), ...
            [mdlName '/Left Heart/Rla'], 'Value', num2str(parameters(44)), ...% heart_params
            [mdlName '/Autonomic Nervous System/Left ventricle/P0lv'], 'Value', num2str(sampleSet(i,22)), ...
            [mdlName '/Autonomic Nervous System/Left ventricle/kElv'], 'Value', num2str(sampleSet(i,23)), ...
            [mdlName '/Autonomic Nervous System/Left ventricle/Vulv'], 'Value', num2str(sampleSet(i,24)), ...
            [mdlName '/Left Heart/Subsystem2/kRlv'], 'Value', num2str(sampleSet(i,25)), ...
            [mdlName '/Autonomic Nervous System/Emaxlv_inv/Emaxlv0_denervated (Emaxlv0)'], 'Value', num2str(sampleSet(i,26)), ...
            [mdlName '/Autonomic Nervous System/Symp. Efferent Pathways/fes_inf'], 'Value', num2str(sampleSet(i,27)), ...
            [mdlName '/Autonomic Nervous System/Carotid Sinus (baroreceptors)/ka'], 'Value', num2str(sampleSet(i,28)), ...
            [mdlName '/Autonomic Nervous System/BP'], 'Value', num2str(BPvalues(j)), ...
            [mdlName '/Autonomic Nervous System/BPswitch'], 'Value', num2str(BPswitch));
    end
end


%% Run simulations
simout = parsim(simIn);
% simtest = sim(simIn(1,1));

%% Analysis
%start & end times for calculating steady state values
tshort = [150 165];     % chose 15 sec range since this is what is used clinically to determine BPM


% pre allocate
simRes = zeros(n,m);
MAP = zeros(n,m);
SP = zeros(n,m);
DP = zeros(n,m);
CO = zeros(n,m);
ECSP = zeros(n,m);
NAactivity = zeros(1,n);
DMVactivity = zeros(1,n);
sympActivity = zeros(1,n);
SDRRstore = zeros(1,n);
RMSSDstore = zeros(1,n);
pNN50store = zeros(1,n);


% Calculate HR and MAP for each parameter set
for i = 1:n
    for j = 1:m
        HR_vals         = simout(i,j).HR_TS;
        % Calculate percent change in heart rate due to step
        if length(mean(getsampleusingtime(HR_vals,165,180),'Weighting','time')) == 0
            simRes(i,j)      = nan; %baseline
        else

            % Calculate hemodynamic metrics
            physOutputs = PhysOutputs_Gen_TS(simout(i,j),tshort,tshort);
            simRes(i,j) = mean(getsampleusingtime(HR_vals,tshort(1),tshort(2)),'Weighting','time');
            SP(i,j) = physOutputs(1);
            DP(i,j) = physOutputs(2);
            MAP(i,j) = (SP(i,j) + 2*DP(i,j))/3;
            CO(i,j) = physOutputs(4);
            if j == zeroIdx
                NAactivity(1,i) = mean(getsampleusingtime(simout(i,j).fevHR,tshort(1),tshort(2)),'Weighting','time');
                DMVactivity(1,i) = mean(getsampleusingtime(simout(i,j).fevEmax,tshort(1),tshort(2)),'Weighting','time');
                sympActivity(1,i) = mean(getsampleusingtime(simout(i,j).fesh,tshort(1),tshort(2)),'Weighting','time');
                % Calculate HRV metrics
                trange = [120 300];
                [SDRR, RMSSD, pNN50] = HRVmetrics_timeDomain(simout(i,j),trange);
                SDRRstore(1,i) = SDRR;
                RMSSDstore(1,i) = RMSSD;
                pNN50store(1,i) = pNN50;
            
            end
        end
    end
    % Convert MAP to ECSP
    ECSP(i,:) = MAP(i,:) + BPvalues;
end
HR = simRes;


%%
% load 'ka[15_15]25-Feb-2024 19_32_17.mat'
%% plot baroreflex curve
% addpath('C:\Users\mmgee\Downloads') %
% addpath('C:\Users\mmgee\Box\Michelle-Gee\Research\MI model')
% addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
% load('C:\Users\mmgee\Downloads\all[15_15]01-Mar-2024 17_51_46.mat')
lw = 2;
ms = 6;
fs = 16;
gray = [0.7 0.7 0.7 0.5];

% load patient data and plot
patientData = {'P1','P2','P3','P4','P5'};

for i = 1:length(patientData)
    load(patientData{i})
end

figure;
for i = 1:length(patientData)
    data = eval(patientData{i});

    % patient data
    x_raw = data.raw_ECSP;
    y_raw = 60./data.raw_HR;
    plot(x_raw,y_raw,'ko','MarkerSize',ms,'MarkerFaceColor','r','HandleVisibility','off')
    x_fit = data.fit_ECSP;
    y_fit = 60./data.fit_HR;
    if i == 1
        plot(x_fit,y_fit,'--','LineWidth',lw,'Color','r');
    else
        plot(x_fit,y_fit,'--','LineWidth',lw,'Color','r','HandleVisibility','off');
    end
    hold on
end

% preallocate vectors to store parameter fits
A1Store = zeros(n,1);
A2Store = zeros(n,1);
A3Store = zeros(n,1);
A4Store = zeros(n,1);

% Check for NaN values and remove rows
nanRowsHR = any(isnan(HR), 2);
nanHRidx = find(nanRowsHR == 1);
nanRowsECSP = any(isnan(ECSP),2);
nanECSPidx = find(nanRowsECSP == 1);
rows2rm = unique([nanECSPidx; nanHRidx]);
mask = true(size(HR, 1), 1);
mask(rows2rm) = false;
HRClean = HR(mask, :);
ECSPClean = ECSP(mask,:);
RR_plot = 60./HRClean;

n_clean = n - length(rows2rm);
% size(ECSPClean)
% ECSPClean
% RR_plot

for i = 1:n_clean
    if i == 1
        plot(ECSPClean(i,:),RR_plot(i,:),'LineWidth',lw/4,'Color',gray)
        hold on
        [A1,A2,A3,A4] = sigmoidRegression(ECSPClean(i,:)',RR_plot(i,:)');
        A1Store(i) = A1;
        A2Store(i) = A2;
        A3Store(i) = A3;
        A4Store(i) = A4;
    else
        plot(ECSPClean(i,:),RR_plot(i,:),'LineWidth',lw/4,'Color',gray,'HandleVisibility','off')
        [A1,A2,A3,A4] = sigmoidRegression(ECSPClean(i,:)',RR_plot(i,:)');
        A1Store(i) = A1;
        A2Store(i) = A2;
        A3Store(i) = A3;
        A4Store(i) = A4;
    end
end



xlabel('Estimated carotid sinus pressure (mm Hg)')
ylabel('RR interval (s)') %
legend('Data from individuals','Individual models','Location','Southeast')
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 600, 500])
hold off
png_filename = [filename '.png'];
saveas(gcf,png_filename)

% reformat sampleSet so it includes all parameters
sampleSetTemp = repmat(parameters,n,1);
sampleSetTemp(:,sigParamCombinedIdx) = sampleSet;
sampleSet = sampleSetTemp;

folder = './';
filename = [folder filename '.mat']
save(filename,'ECSP','HR','sampleSet','SP', 'DP','CO', 'MAP','A1Store','A2Store','A3Store','A4Store','NAactivity','DMVactivity','sympActivity','SDRRstore', 'RMSSDstore', 'pNN50store')
end