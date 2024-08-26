clear; close all;
%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);

%% Add file paths
my_dir = pwd;
addpath(genpath(my_dir))
addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')

%% Main
% control PI parameters
PNDMVdelay = 0.3;
ICNparams = [1.6906 7.3295 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 PNDMVdelay 3.329861 2.661685 5.642977 0.066794];

% NA/DMV
NAparams = [4.88, 15.78, 59.83, 23, 0.61, 11, 12.81, 7, 2.5901, 6.66, 42.91, 33.5]; %  NA, NActr, DMV

% NTS
NTSparams = [0.30, 21.50, 37.07, 21, 0.45, 28.33, 10.2, 7, 2.75, 31.57, 11.13, 2]; %  BR, CPR, LSR; fmin, fmax, fmid, k

% heart_params
Cla = 19.23;
Vula = 25;
Rla = 2.5e-3;
P0lv = 1.5;
kElv = 0.014;
Vulv = 16.77;
kRlv = 3.75e-4;
Emaxlv0 = 1.283;
fes_inf = 2.1;
BPvalues = -70:10:90; % neck chamber pressure


ParkHealthy = [Cla, Vula, Rla,P0lv,kElv,Vulv,kRlv, Emaxlv0, fes_inf];
ka = 11.758;

% variables that change
parameters = [ICNparams, NTSparams, NAparams, ParkHealthy, ka]; % expand this later
filename = 'all[25_25]500_06182024';


[sampleSet, SP, DP, HR, MAP,A1Store,A2Store,A3Store,A4Store,NAactivity,DMVactivity,sympActivity] = barocurve_dist(BPvalues,parameters,filename);


%% find indices of sample sets that are within the range of observed
% baroreflex curves
filename = 'all[25_25]500_06182024';
addpath('C:\Users\mmgee\Downloads') %
addpath('C:\Users\mmgee\Box\Michelle-Gee\Research\MI model')
addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
load([filename '.mat'])
% load('C:\Users\mmgee\Box\Michelle-Gee\Research\MI model\all[15_15]50003-Mar-2024 20_11_04.mat')
% addpath(genpath('C:\Users\mmgee\MATLAB\umapFileExchange (4.4)'))

%% Identify sample sets where there was a nan value for HR/ECSP because A1-A4
% was not calculated
nanRowsHR = any(isnan(HR), 2);
nanHRidx = find(nanRowsHR == 1);
nanRowsECSP = any(isnan(ECSP),2);
nanECSPidx = find(nanRowsECSP == 1);
rows2rm = unique([nanECSPidx; nanHRidx]);
mask = true(size(HR, 1), 1);
mask(rows2rm) = false;

% Remove rows
HRstore = HR(mask, :);
ECSPstore = ECSP(mask,:);
COstore = CO(mask,:);
DMVactivityStore = DMVactivity(:,mask)';
DPstore = DP(mask,:);
MAPstore = MAP(mask,:);
NAactivityStore = NAactivity(:,mask)';
RMSSD = RMSSDstore(:,mask)';
SDRR = SDRRstore(:,mask)';
SPstore = SP(mask,:);
PNN50 = pNN50store(:,mask)';
sampleSetStore= sampleSet(mask,:);
sympActivityStore = sympActivity(:,mask)';


% Remove sample sets where A1-A4 are zero (due to nan HR value)
A = [A1Store A2Store A3Store A4Store];
[rowsRm, colsRm] = find(A == 0);
mask = true(size(A, 1), 1);
mask(rowsRm) = false;
A = A(mask, :);

% Unpack A1-4 back into column vectors
A1Store = A(:,1);
A2Store = A(:,2);
A3Store = A(:,3);
A4Store = A(:,4);

disp('size of A')
disp(size(A))

disp('size of sample set all')
disp(size(sampleSetStore))

matFilename = 'FilteredDistribution_061824.mat';
save(matFilename, 'A1Store', 'A2Store', 'A3Store', 'A4Store', 'COstore','DMVactivityStore', ...
    'DPstore', 'ECSPstore', 'HRstore','MAPstore','NAactivityStore','RMSSD','PNN50', ...
    'SDRR', 'SPstore','sampleSetStore','sympActivityStore')

%% find indices of sample sets that are within the range of observed
load(matFilename)
BPvalues = -70:10:90; % neck chamber pressure
zeroIdx = find(BPvalues == 0);

% 95% confidence intervals for these parameters from Seredynski. Note that
% the range for A2 is bounded by zero for the higher value so that the
% slope of the sigmoid is not negative
A1range = [0.0316    0.4639];
A2range = [-0.3006    0];
A3range = [60.0981  128.6592];
A4range = [0.5480    1.0360];
RRrange = [0.624779 1.16102]; % min and max reported RR interval values
RR_plot = 60./HRstore;



[rows, cols] = size(sampleSetStore);
idxStore = [];
% mastitskaya HR range 67.2-76.8, MAP 59.8-122.2
count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
count5 = 0;

% A1 positive, A2 negative
% Use baseline average HR and MAP from Seredynski and Mastitskaya to
% convert between human and rat
% Use larger standard deviations (Mastitskaya) to determine ranges
% calculated in Mastitskaya_95CI.m
HRmean = 65.4448;
MAPmean = 99.5180;
count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
count5 = 0;
count6 = 0;
for i = 1:rows
    if (A1Store(i) > 0 && A2Store(i) < 0) || (A1Store(i) < 0 && A2Store(i) > 0)
        count1 = count1 + 1;
        if abs(A1Store(i)) > A1range(1) && abs(A1Store(i)) < A1range(2)
            count2 = count2 + 1;
            if -abs(A2Store(i)) > A2range(1) && -abs(A2Store(i)) < A2range(2)
                count3 = count3 + 1;
                if A3Store(i) > A3range(1) && A3Store(i) < A3range(2)
                    count4 = count4 + 1;
                    if A4Store(i) > A4range(1) && A4Store(i) < A4range(2)
                        count5 = count5 + 1;
                        if HRstore(i,zeroIdx) > 55 && HRstore(i,zeroIdx) < 76 % filter baseline HR based on 95% CI from Mastitskaya
                            count6 = count6 + 1;
                            if MAPstore(i,zeroIdx) > 69 && MAPstore(i,zeroIdx) < 130 % filter baseline MAP based on 95% CI from Mastitskaya
                                % if RR_plot(i,1) > 0.6 && RR_plot(i,end) < 1.2
                                    idxStore = [idxStore i];
                                % end
                            end
                        end
                    end
                end
            end
        end
    end
end

count1
count2
count3
count4
count5
count6


% idxStore([5,41]) = [];
save('FilteredDistributionIdx_061824.mat','idxStore','sampleSet')
%% load patient data
patientData = {'P1','P2','P3','P4','P5'};

for i = 1:length(patientData)
    load(patientData{i})
end

lw = 2;
ms = 6;
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

% ECSP = MAP(:,zeroIdx) + BPvalues;
gray = [0.7 0.7 0.7 0.5];
fs = 14;


for i = 1:length(idxStore)
    if i == 1
        plot(ECSP(idxStore(i),:),RR_plot(idxStore(i),:),'LineWidth',lw/4,'Color',gray)
    else
        plot(ECSP(idxStore(i),:),RR_plot(idxStore(i),:),'LineWidth',lw/4,'Color',gray,'HandleVisibility','off')
    end
end



xlabel('Estimated carotid sinus pressure (mm Hg)')
ylabel('RR interval (s)') %
legend('Data from individuals','Individual models','Location','Southeast')
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 600, 500])
png_filename = 'FilteredDistribution061824.png';
saveas(gcf,png_filename)
% matFilename = [filename 'FilteredDistribution061824.mat'];
% save(matFilename,'idxStore','sampleSet','MAP','HR')
disp('Number of accepted parameter sets')
disp(length(idxStore))

%% Clustering of accepted preIR model properties

% deep breath
breathLength = 5.8;
parameters = sampleSet(idxStore,:);
[RSAamp, stdRSA] = deepBreath(breathLength,parameters);

% HRV
% [SDRR, RMSSD, LFHFratio,pRMS,powbp_total,powbp_vlf,powbp_lf,powbp_hf] = HRVmetrics(simOut,trange);

% baseline hemodynamic variables (from .mat)
%% Calculate NA and DMV activity for 29 accepted preIR models
% load 'FilteredDistribution031924.mat' %sampleSet, idxStore
%
% parameters = sampleSet(idxStore,:);
% filename = 'preIRvagalActivity_031924';
% [NAactivity,DMVactivity,sympActivity] = calcVagalActivity(parameters,filename);
%% tSNE of parameters
% z score data
% Z = zscore(sampleSet); % columns centered to zero mean and SD is 1
% 
% numParamSets = rows;
% labels = repmat({'Rejected'},numParamSets,1); % (223,numParamSets,1); %
% labels(idxStore) = {'Accepted'};
% % umap_input = [Z(:,[18:34]) labels];
% 
% rng default
% perplexity = 35;
% % end_col = size(umap_input,2);
% Y = tsne(sampleSet,'Algorithm','exact','Distance','cityblock','Perplexity',perplexity,'Standardize',1);
% % [reduction, umap, clusterIdentifiers, extras]=run_umap(umap_input,'marker_size',10,'save_output','true','label_column',end_col)
% 
% figure;
% 
% gscatter(Y(:,1),Y(:,2),labels,[0.7 0.7 0.7;0 0 1]);
% xlabel('t-SNE 1')
% ylabel('t-SNE 2')
% title(num2str(perplexity))
% filename = 'preIR_params_tSNE_04152024.png'; %
% set(gca,'FontSize',16)
% saveas(gcf,filename)


%% Pair plots
% acceptedIndices = idxStore;
% varNames = {'PN_{NA} f_{min}','PN_{NA} f_{max}','PN_{NA} f_{mid}','PN_{NA} k', ...
%     'LCN f_{min}','LCN f_{max}','LCN f_{mid}','LCN k', ...
%     'PN_{DMV} f_{min}','PN_{DMV} f_{max}','PN_{DMV} f_{mid}','PN_{DMV} k', '\tau_{PNDMV}', ...
%     'k_{fevEmax}','k_{BR}','k_{CP}','k_{fesh}', ...
%     'BR f_{min}','BR f_{max}','BR f_{mid}','BR k', ...
%     'CP f_{min}','CP f_{max}','CP f_{mid}','CP k', ...
%     'LSR f_{min}','LSR f_{max}','LSR f_{mid}','LSR k', ...
%     'NA f_{min}','NA f_{max}','NA f_{mid}','NA k', ...
%     'NActr f_{min}','NActr f_{max}','NActr f_{mid}','NActr k', ...
%     'DMV f_{min}','DMV f_{max}','DMV f_{mid}','DMV k', ...
%     'C_{la}','P_{0lv}','k_{Elv}','V_{ulv}','fes_{inf}','k_a'};
% filename = 'all[2_2]_pairplot.png';
% generatePairPlots(sampleSet, acceptedIndices, varNames, filename)
%% Exit code
close_system
delete(myPool);
exit
