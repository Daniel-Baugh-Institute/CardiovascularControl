clear; close all;
%%
addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')
cd ..
my_dir = pwd;
addpath(genpath(my_dir))
%% Set up parallel pool
% myCluster = parcluster('local');
% myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
% myCluster.JobStorageLocation = getenv('TMPDIR');
% myPool = parpool(myCluster, myCluster.NumWorkers);

%% Fit Osculati 1990 mean data
% load data
neckChamberPressure = [8	13.2	18.5	23.7	30	37.5]; % mm hg
RRbase = 777.8; % +/- 34.7 msec
deltaRR = [46.50360279	64.79801636	89.04640376	126.5452133	144.0065953	171.9898378];
SEM = [7.550265772	10.06205045	13.11093148	24.29514847	24.91603345	24.7516633]; % msec
num = 12;
SD = SEM.* sqrt(num); % standard deviation
CI95 = 1.96.* SD;
RR_Osculati = RRbase + deltaRR;

% Curve not a good fit for data
% [A1,A2,A3,A4] = sigmoidRegression(neckChamberPressure',RR') % inputs must be column vectors
%% Plot
% data
x = neckChamberPressure;
y = RR_Osculati;
errors = CI95;

% Plot the line
ms = 8;
figure;
hold on;
errorbar(x, y, errors,'ro', 'MarkerFaceColor','r','MarkerSize',ms);

% Add labels and title
xlabel('Neck chamber pressure (mm Hg)');
ylabel('RR interval (ms)');
% title('Line with Shaded Region (Standard Error)');
legend('Osculati 1990')
set(gca,'FontSize',14)

% Show plot
hold off;

%% change these liens
% load PI individual parameters
% matFilename = [filename 'FilteredDistribution.mat'];
date = '06232024';
load('FilteredDistribution_061824.mat')
matFilename = 'FilteredDistributionIdx_061824.mat';
load(matFilename) % gives 'idxStore','sampleSet'


HR_preIR = HRstore;
MAP = MAPstore;
NAactivity_preIR = NAactivityStore;
DMVactivity_preIR = DMVactivityStore;
sympActivity_preIR = sympActivityStore;



% Change these
numPatients = length(idxStore);


%% change this line
% countTotal = 0;
% BPvalues = -70:10:90; % neck chamber pressure preIR
% zeroIdx = find(BPvalues == 0);
% BPvalues_postIR = [-70 -50 -25 0 8	13.2	18.5	23.7	30	37.5 50 75 90]; 
% 
% mdlName     = 'ICN_model_v15_Mastitskaya2012_control_r2020b';
% % idxFail = 212;%[1 2 3 4 5 6 7 73 100 107 134 138 150 167];
% for i = 36:40%numPatients
% 
%     parameters = sampleSet(idxStore(i),:);
%     numSampleSets = 100;
%     variationFactor = 1.5;
%     mdlAlternative = 6;
%     filename = ['postIR_all_' date '_' num2str(i)];
%     HR_PI_mdl = HR_preIR(idxStore(i),zeroIdx);
%     MAP_PI_mdl = MAP(idxStore(i),zeroIdx);
% 
%     BRS_preIR = A2Store(idxStore(i)); % preIR sigmoid slope parameter
%         [HR_mdl,MAP_mdl,SP_postIR,DP_postIR,CO_postIR,paramMatrix,...
%     NA_activity_postIR,DMV_activity_postIR, sympActivity, SDRR_postIR, RMSSD_postIR, ...
%     pNN50_postIR, crit, A] = paramEstimationHRMAP_global(parameters,...
%     mdlName,numSampleSets,variationFactor,HR_PI_mdl,MAP_PI_mdl,...
%     mdlAlternative,NAactivity_preIR,DMVactivity_preIR,...
%     sympActivity_preIR,BRS_preIR,filename);
% 
% 
%     save(['./plots/paramEstimation/' filename '.mat'], ...
%         'HR_mdl','MAP_mdl','SP_postIR','DP_postIR','CO_postIR','paramMatrix',...
%     'NA_activity_postIR','DMV_activity_postIR', 'sympActivity', 'SDRR_postIR', 'RMSSD_postIR', ...
%     'pNN50_postIR', 'crit', 'A')
%     i
% 
%     % Combine all data into a struct for all preIR models
%     postIRstore(i).MAP = MAP_mdl; 
%     postIRstore(i).HR = HR_mdl;
%     postIRstore(i).SP = SP_postIR;
%     postIRstore(i).DP = DP_postIR;
%     postIRstore(i).CO = CO_postIR;
%     postIRstore(i).HR = HR_mdl;
%     postIRstore(i).params = paramMatrix;
%     postIRstore(i).NA = NA_activity_postIR;
%     postIRstore(i).DMV = DMV_activity_postIR;
%     postIRstore(i).symp = sympActivity;
%     postIRstore(i).SDRR = SDRR_postIR;
%     postIRstore(i).RMSSD = RMSSD_postIR;
%     postIRstore(i).pNN50 = pNN50_postIR;
%     postIRstore(i).crit = crit;
%     postIRstore(i).Avals = A;
% end
%%
matFilename = ['postIR_all_' date '.mat'];
load(matFilename)
% save(matFilename, 'postIRstore') % postIR_idx is model sub idx for accepted postIR models

%% Unpack postIR models to plot HRMAP bar plot for those that passed HR, MAP criteria
% preIR_MAP_models = [];
% preIR_HR_models = [];
% postIR_MAP_models = [];
% postIR_HR_models =[];
% 
% % indices of non-empty fields
% indices_i = [];
% indices_j = [];
% postIR_HR_models = [];
% postIR_MAP_models = [];
% preIR_idx = [];
% reps = zeros(1,numPatients);
% 
% % Loop through struct array
% for i = 1:numPatients % preIR idx loop
%     indices_j = [];
%     for j = 1:numSampleSets % postIR idx loop
%     % Check if systolic pressure is in range
%         if postIRstore(i).crit(j,1) == 1
%             % Check HR
%             if postIRstore(i).crit(j,2) == 1
%                 % Check MAP decrease
%                 % if postIRstore(i).crit(j,3) == 1
%                     % Check HR within 20% change
%                     % if postIRstore(i).crit(j,4) == 1
%                         indices_i = [indices_i, i];
%                         indices_j = [indices_j, j]; 
%                         postIR_HR_models = [postIR_HR_models postIRstore(i).HR(j)];
%                         postIR_MAP_models = [postIR_MAP_models postIRstore(i).MAP(j)];
%                     % end
%                 % end
%             end
%         end
%     end
%     reps(i) = length(indices_j); % number of accepted postIR mdoels per preIR model
%     preIR_idx = [preIR_idx i];
% end
% 
% 
% for i = 1:numPatients
%     preIR_HR_models = [preIR_HR_models repmat(HR_preIR(idxStore(preIR_idx(i)),zeroIdx),1,reps(i))];
%     preIR_MAP_models = [preIR_MAP_models repmat(MAP(idxStore(preIR_idx(i)),zeroIdx),1,reps(i))];
% end
% 
% altMdlName = 'all';
% Mastitskaya_plot_bar_HRMAP(preIR_MAP_models, postIR_MAP_models, preIR_HR_models, postIR_HR_models,altMdlName)

%% Plot parameter distributions for pre-IR and post-IR
% preIR_parameters = sampleSet(idxStore,42:46);
% postIR_parameters = paramStore(:,42:46);
% parameterNames = {'C_{la} (mL/mm Hg)','P_{0lv} (mm Hg)','k_{Elv} (mm Hg^{-1})','V_{ulv} (mL)','fes_{inf} (Hz)'};
% saveFilename = ['all_paramDist_' date '.png'];
% plot_paramDist(preIR_parameters, postIR_parameters, parameterNames, saveFilename)

%% tSNE of crit matrix
% % color by preIR model?-- yes to tell us which individuals can and can't
% % compensate
% perplexity = 5;
% rng default
% 
% % format data into a single matrix
% Xdata = zeros(numPatients*numSampleSets,7);
% for i = 1:numPatients
%     start = 1 + (i-1)*numSampleSets;
%     stop = i*numSampleSets;
% 
%     Xdata(start:stop,1:7) = postIRstore(i).crit(:,1:7);
% end
% 
% 
% % tSNE-- label postIR vs preIR
% Y = tsne(Xdata,'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',0);% 
% 
% % label based on preIR model
% pre_idx = 1:1:size(Xdata,1);
% % preallocate
% % obslabels = repmat({'No post-IR model'},1,size(Xdata,1));
% % obslabels(hasPostIRmdlIdx) = {'Has post-IR model'};
% obslabels = repmat({num2str(1)},numPatients*numSampleSets,1);
% 
% for i = 2:numPatients
%     start = 1 + (i-1)*numSampleSets;
%     stop = i*numSampleSets;
%     labels = {num2str(i)};
%     idx = start:1:stop;
%     obslabels(idx) = labels;
% end
% 
% 
% 
% figure;
% b = gscatter(Y(:,1),Y(:,2),obslabels);
% xlabel('t-SNE 1')
% ylabel('t-SNE 2')
% set(gca,'FontSize',16)
% saveas(gcf,'prepostIRaccepted_tsne.png')
%% Plot accepted baroreflex curves

postIRstore_filename = 'postIR_all_06232024';
postIRbarocurve_filename = 'postIR_all_06182024';
BPvalues = -70:10:90;
brown = [196/255 153/255 124/255];
color = brown;
filename = ['overlaid_barocurve_all' date '.png'];
plot_overlaid_barocurve(postIRstore_filename,color,postIRbarocurve_filename)


%% Exit code
% close_system
% delete(myPool);
% exit