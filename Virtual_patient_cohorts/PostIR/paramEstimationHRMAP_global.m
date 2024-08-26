function [HR_mdl,MAP_mdl,SP_postIR,DP_postIR,CO_postIR,paramMatrix,...
    NA_activity_postIR,DMV_activity_postIR,sympActivity,SDRR_postIR, RMSSD_postIR, ...
    pNN50_postIR, crit, A] = paramEstimationHRMAP_global(parameters,...
    mdlName,numSampleSets,variationFactor,HR_PI_mdl,MAP_PI_mdl,...
    mdlAlternative,NA_activity_preIR,DMV_activity_preIR,...
    sympActivity_preIR,BRS_preIR,filename)
% Tuning HR and MAP to fit Mastitskaya 2012 data for any alternative model
% (cardiac, baroreceptors,NTS, NA/DMV, ICN)
% mdlAlternative = 1, cardiac
% mdlAlternative = 2, baroreceptors
% mdlAlternative = 3, NTS
% mdlAlternative = 4, NA/DMV
% mdlAlternative = 5, ICN
% target: [HR; MAP] target values

% Michelle Gee
% March 1, 2024

%% CHANGES:
% save all data for all models
% output crit 1-7 as matrix
% make sure I'm okay with filtering criteria (using Mastitskaya as
% qualitative data
% Save A values in matrix

%%
addpath('C:\Users\mmgee\MATLAB\Projects\untitled')
%% Set up parallel pool
% myCluster = parcluster('local');
% myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
% myCluster.JobStorageLocation = getenv('TMPDIR');
% myPool = parpool(myCluster, myCluster.NumWorkers);


%% Create parameter set and number of simulations


% Parameters that vary
paramVec = parameters;
BPswitch = -1; % closed loop
kRSA = 0.5;

% Set parameter range for Sobol sampling based on alternative models
factor = variationFactor;
if mdlAlternative == 1
    paramVec = parameters([42,45:50]);
    minParams = paramVec./factor;
    maxParams = paramVec.*factor;
elseif mdlAlternative == 2 % baroreceptor model
    paramVec = parameters([42,45:51]);
    minParams = paramVec./factor; % bounded ka so baroreflex gain lower
    minParams(6) = parameters(51);
    maxParams = paramVec.*factor;
elseif mdlAlternative == 3 % NTS model
    paramVec = [parameters([20 23 27 28]) parameters([42,45:50])]; % CP_fmax (23), BR_fmid (20), LS_fmid (28), LS_fmax (27)
    minParams = paramVec./factor;
    % restrict to increased NTS output (fmax)
    minParams(2) = paramVec(2);
    minParams(3) = paramVec(3);

    maxParams = paramVec.*factor;
    % restrict to lower fmid
    maxParams(1) = paramVec(1);
    maxParams(4) = paramVec(4);
elseif mdlAlternative == 4 % NADMV
    paramVec = [parameters([30, 31, 32, 36, 39, 40]) parameters([42, 45:50])]; % NA_fmin (30), NA_fmid (32), DMV_fmid (40), DMV_fmax (39), NActr_fmid (36)
    minParams = paramVec./factor;
    maxParams = paramVec.*factor;
elseif mdlAlternative == 5 % ICN
    paramVec = [parameters([1,2,3,4,9,11]) parameters([42,45:50])]; % PNDMV_fmin (9), PNNA_fmid (3), PNNA_fmin (1), PNNA_fmax (2), PNDMV_fmid (11), PNNA_k (4)
    minParams = paramVec./factor;
    maxParams = paramVec.*factor;
elseif mdlAlternative == 6
    paramVec = parameters([1 2 3 4 5 7 8 9 11 16 20 23 27 28 30 31 32 36 39 40 42 45 46 47 48 49 50 51]);
    minParams = paramVec./factor;
    maxParams = paramVec.*factor;

    % additional constraints for BR and NTS
    minParams(end) = parameters(51);
    maxParams(11) = paramVec(11);
    minParams(12) = paramVec(12);
    maxParams(13) = paramVec(13);
    minParams(14) = paramVec(14);

else
    disp('Not a valid model alternative number')
    return
end

[sobolSet] = SobolSample(minParams,maxParams,length(paramVec),numSampleSets);

% insert sobol sampled values into parameter matrix based on alternative
% model
sampleSet = repmat(parameters,numSampleSets,1);
if mdlAlternative == 1 % cardiac
    sampleSet(:,[42,45:50]) = sobolSet;
elseif mdlAlternative == 2 % baroreceptors
    sampleSet(:,[42,45:51]) = sobolSet;
elseif mdlAlternative == 3 % NTS
    sampleSet(:,[20 23 27 28]) = sobolSet(:,1:4);
    sampleSet(:,[42,45:50]) = sobolSet(:,5:11);
elseif mdlAlternative == 4 % NA/DMV
    sampleSet(:,[30, 31, 32, 36, 39, 40]) = sobolSet(:,1:6);
    sampleSet(:,[42,45:50]) = sobolSet(:,7:13);
elseif mdlAlternative == 5 % ICN
    sampleSet(:,[1,2,3,4,9,11]) = sobolSet(:,1:6);
    sampleSet(:,[42,45:50]) = sobolSet(:,7:13);
else % all parameters
    sampleSet(:,[1 2 3 4 5 7 8 9 11 16 20 23 27 28 30 31 32 36 39 40 42 45 46 47 48 49 50 51]) = sobolSet;
end
sampleSet_postIR = sampleSet;

% BPvalues = -70:16:90;% from main_ACC: BPvalues = -70:5:90;
% BPvals = BPvalues';
% m = length(BPvals);






%% plot baroreflex curve and check criteria
BPvalues = [-70 -50 -25 0 8	13.2	18.5	23.7	30	37.5 50 75 90]; % mm hg. must include 0 for function
% "negative neck chamber pressure effectively increases carotid sinus
% pressure" Osculati


folder = './plots/plots20240410/';

A1store = zeros(numSampleSets,1);
A2store = zeros(numSampleSets,1);
A3store = zeros(numSampleSets,1);
A4store = zeros(numSampleSets,1);

filenamePNG = [folder filename '_barocurve'];

% control PI parameters
ICNparams = sampleSet(:,1:17);%)[0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 PNDMVdelay 3.329861 2.661685 5.642977 0.066794];

% NA/DMV
NAparams = sampleSet(:,30:41);%[4.88, 15.78, 59.83, 23, 0.61, 11, 12.81, 7, 2.5901, 6.66, 42.91, 33.5]; %  BR, CPR, LSR

% NTS
NTSparams = sampleSet(:,18:29);%[0.30, 21.50, 37.07, 21, 0.45, 28.33, 10.2, 7, 2.75, 31.57, 11.13, 2]; %  BR, CPR, LSR; fmin, fmax, fmid, k

params = sampleSet(:,[42,45:51]);
heart_params = params(:,1:7);
size(heart_params)
ka = params(:,8);
disp('ka length')
length(ka)

baseline = 1; % calculate baseline hemodynamic measurements
[matFilename, HRvalidation, MAP_EI_mdl, SP_store, DP_store, CO_store, ...
    NA_activity_postIR, DMV_activity_postIR, symp_activity_postIR, ...
    SDRR_store, RMSSD_store, pNN50_store] = baroreflex_curve(mdlName, ...
    BPvalues,filename,ICNparams,NAparams,NTSparams,heart_params,ka,baseline);

% returns matfile with variables 'SBP','deltaRR','BPvalues','HR','baselineMAP'
load(matFilename,'SBP','deltaRR','BPvalues','HR','baselineMAP')


%% Fit baroreflex curve
for i = 1:numSampleSets
    ECSP = baselineMAP + BPvalues;
    RR = 60./HR(i,:)';
    if any(isnan(ECSP)) || any(isnan(RR))
        idx2rmECSP = find(isnan(ECSP));
        idx2rmRR = find(isnan(RR));
        idx2rm = union(idx2rmRR,idx2rmECSP);
        RRclean = RR;
        RRclean(idx2rm) = [];
        ECSPclean = ECSP;
        ECSPclean(idx2rm) = [];
        if length(ECSPclean) <= 4
            A1 = nan;
            A2 = nan;
            A3 = nan;
            A4 = nan;
            disp('Warning: less than 4 data points in baroreflex curve simulation so sigmoidal model was not fit')
        else
            [A1,A2,A3,A4] = sigmoidRegression(ECSPclean',RRclean);
        end
    else
        [A1,A2,A3,A4] = sigmoidRegression(ECSP',RR);
    end
    A1store(i) = A1;
    A2store(i) = A2;
    A3store(i) = A3;
    A4store(i) = A4;
end

% Combine A1-A4
A = [A1store A2store A3store A4store];
%% Accept or reject parameter sets based on the following criteria:
% 1. within the experimentally observed MAP range
% 2. within the experimentally observed HR range
% 3. decrease in MAP from preIR to postIR
% 4. less than 3% change in HR
% 5. within the range observed in the baroreflex curve (do this after
% filtering the first two criteria so we don't have to run the model as
% many times)


% Osculati data for comparison
RRbase = 777.8; % +/- 34.7 msec
deltaRR = [46.50360279	64.79801636	89.04640376	126.5452133	144.0065953	171.9898378];
SEM = [7.550265772	10.06205045	13.11093148	24.29514847	24.91603345	24.7516633]; % msec
num = 12;
SD = SEM.* sqrt(num); % standard deviation
CI95 = 1.96.* SD;
RR_Osculati = RRbase + deltaRR;
minOsculati = RR_Osculati - CI95;
maxOsculati = RR_Osculati + CI95;


% postIR model subset index, hemodynamic constraints, vagal activity decrease, BR slope decrease
crit1 = zeros(numSampleSets,1);
crit2 = zeros(numSampleSets,1);
crit3 = zeros(numSampleSets,1);
crit4 = zeros(numSampleSets,1);
crit5 = zeros(numSampleSets,1);
crit6 = zeros(numSampleSets,1);
crit7 = zeros(numSampleSets,1);
all_crit = zeros(numSampleSets,1);

% CHANGED HERE
symp_vagal_ratio_preIR = sympActivity_preIR/(NA_activity_preIR + DMV_activity_preIR);


for i = 1:numSampleSets % 10
    % Changed to test if ratio of parasymp/symp activity decrease vs
    % decrease in parasymp activity criteria makes a difference
    % NAchange = (NA_activity_postIR(i) - NA_activity_preIR)/NA_activity_preIR;
    % DMVchange = (DMV_activity_postIR(i) - DMV_activity_preIR)/DMV_activity_preIR;


    % Systolic pressure in range of Osculati
    % Osculati SBP mean +/- SD = 118.7750 +/- 11.6290; 95% CI: [95.99, 141.57]
    if SP_store(i) > 95.99 & SP_store(i) < 141.57%MAP_EI_mdl(i) > 62 & MAP_EI_mdl(i) < 98 %MAP_EI_mdl(i) > 60.3 & MAP_EI_mdl(i) < 95.9
        crit1(i) = 1;
    end

    % Heart rate in range of Osculati
    if HRvalidation(i) > 58 & HRvalidation(i) < 116 %HRvalidation(i) > 54 & HRvalidation(i) < 70 %HRvalidation(i) > 61.5 & HRvalidation(i) < 75.5
        crit2(i) = 1;
    end

    % MAP decrease based on Mastitskaya
    if MAP_EI_mdl(i) <= MAP_PI_mdl
        crit3(i) = 1;
    end

    % HR relatively unchanged based on Mastitskaya
    if abs(HRvalidation(i) - HR_PI_mdl)/HR_PI_mdl < 0.2
        crit4(i) = 1;
    end

    % Sympathetic to vagal activity ratio
    symp_vagal_ratio_postIR = symp_activity_postIR(i)/(NA_activity_postIR(i) + DMV_activity_postIR(i));
    if symp_vagal_ratio_postIR > symp_vagal_ratio_preIR %NAchange < 0.1
        crit5(i) = 1;
    end

    % Check if baroreflex curve is within Osculati 95% confidence interval
    RRmdl = 1000.*60./HR(i,:);
    withinRange = RRmdl(5:10) >= minOsculati' & RRmdl(5:10) <= maxOsculati';
    if sum(withinRange) == 6
        crit6(i) = 1;
    end

    % Check if sigmoid slope (A2) decreased
    if abs(A2store(i)) < abs(BRS_preIR)
        crit7(i) = 1;
    end

    % Sum models that pass all criteria
    all_crit_sum = crit1 + crit2 + crit3 + crit4 + crit5 + crit6 + crit7;
    if all_crit_sum == 7
        all_crit(i) = 1;
    end

end

%Combine all criteria
crit = [crit1 crit2 crit3 crit4 crit5 crit6 crit7 all_crit];
% Print numbers of models that passed criteria
disp('Number of models that passed all criteria: ')
disp(sum(all_crit))

disp('Criteria 1 - criteria 6')
sum(crit1)
sum(crit2)
sum(crit3)
sum(crit4)
sum(crit5)
sum(crit6)

% Format function outputs
HR_mdl = HRvalidation;
MAP_mdl = MAP_EI_mdl;
SDRR_postIR = SDRR_store;
RMSSD_postIR= RMSSD_store;
pNN50_postIR = pNN50_store;
SP_postIR = SP_store;
DP_postIR = DP_store;
CO_postIR = CO_store;
paramMatrix = sampleSet_postIR;
sympActivity = symp_activity_postIR;

end  