%% Compare the pre- and post-IR relative influence of central (NTS, NA, DMV) and peripheral (ICN, baroreceptors) on heart rate
% Approach: Regression where the inputs are the parameters and the output
% is heart rate. Average the pre-IR central parameter Shapley values and
% compare to the pre-IR peripheral Shapley values. Then compare how this
% ratio changes pre- to post-IR
% clear;
close all;
%% Load data
% addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
% addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
% addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')
my_dir = pwd;
addpath(genpath(my_dir))

BPvalues = -70:10:90; % neck chamber pressure preIR
zeroIdx = find(BPvalues == 0);

% preIR data
load 'FilteredDistribution_061824.mat'
paramsTemp_preIR = sampleSetStore;
HR_preIR = HRstore;

load 'FilteredDistributionIdx_061824.mat'
numPatients = length(idxStore);
numSamples = 100;

% load postIR data
paramsTemp_BR = [];
HRstore_BR = [];
COstore_BR = [];
crit_BR = [];
A2_BR = [];
NA_BR = [];
DMV_BR = [];
symp_BR = [];
preIRmdl_BR = [];
BR_idx=[];

paramsTemp_ICN = [];
HRstore_ICN = [];
COstore_ICN = [];
crit_ICN = [];
A2_ICN = [];
NA_ICN = [];
DMV_ICN = [];
symp_ICN = [];
preIRmdl_ICN = [];
ICN_idx=[];

paramsTemp_cardiac = [];
HRstore_cardiac = [];
COstore_cardiac = [];
crit_cardiac = [];
A2_cardiac = [];
NA_cardiac = [];
DMV_cardiac = [];
symp_cardiac = [];
preIRmdl_cardiac = [];
cardiac_idx=[];

paramsTemp_NADMV = [];
HRstore_NADMV = [];
COstore_NADMV = [];
crit_NADMV = [];
A2_NADMV = [];
NA_NADMV = [];
DMV_NADMV = [];
symp_NADMV = [];
preIRmdl_NADMV = [];
NADMV_idx=[];

paramsTemp_NTS = [];
HRstore_NTS = [];
COstore_NTS = [];
crit_NTS = [];
A2_NTS = [];
NA_NTS = [];
DMV_NTS = [];
symp_NTS = [];
preIRmdl_NTS = [];
NTS_idx=[];

paramsTemp_all = [];
HRstore_all = [];
COstore_all = [];
crit_all = [];
A2_all = [];
NA_all = [];
DMV_all = [];
symp_all = [];
preIRmdl_all = [];
all_idx=[];


for i = 1:numPatients
    for j = 1:numSamples
        load 'postIR_baroreceptors_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_BR = [paramsTemp_BR; postIRstore(i).params(j,:)];
            HRstore_BR = [HRstore_BR postIRstore(i).HR(j)];
            COstore_BR = [COstore_BR postIRstore(i).CO(j)];
            crit_BR = [crit_BR postIRstore(i).crit(j,6)];
            A2_BR = [A2_BR postIRstore(i).Avals(j,2)];
            NA_BR = [NA_BR postIRstore(i).NA(j)];
            DMV_BR = [DMV_BR postIRstore(i).DMV(j)];
            symp_BR = [symp_BR postIRstore(i).symp(j)];
            preIRmdl_BR = [preIRmdl_BR i];
            BR_idx=[BR_idx j];
        end


        load 'postIR_cardiac_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_cardiac = [paramsTemp_cardiac; postIRstore(i).params(j,:)];
            HRstore_cardiac = [HRstore_cardiac postIRstore(i).HR(j)];
            COstore_cardiac = [COstore_cardiac postIRstore(i).CO(j)];
            crit_cardiac = [crit_cardiac postIRstore(i).crit(j,6)];
            A2_cardiac = [A2_cardiac postIRstore(i).Avals(j,2)];
            NA_cardiac = [NA_cardiac postIRstore(i).NA(j)];
            DMV_cardiac = [DMV_cardiac postIRstore(i).DMV(j)];
            symp_cardiac = [symp_cardiac postIRstore(i).symp(j)];
            preIRmdl_cardiac = [preIRmdl_cardiac i];
            cardiac_idx=[cardiac_idx j];
        end


        load 'postIR_NADMV_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_NADMV = [paramsTemp_NADMV; postIRstore(i).params(j,:)];
            HRstore_NADMV = [HRstore_NADMV postIRstore(i).HR(j)];
            COstore_NADMV = [COstore_NADMV postIRstore(i).CO(j)];
            crit_NADMV = [crit_NADMV postIRstore(i).crit(j,6)];
            A2_NADMV = [A2_NADMV postIRstore(i).Avals(j,2)];
            NA_NADMV = [NA_NADMV postIRstore(i).NA(j)];
            DMV_NADMV = [DMV_NADMV postIRstore(i).DMV(j)];
            symp_NADMV = [symp_NADMV postIRstore(i).symp(j)];
            preIRmdl_NADMV = [preIRmdl_NADMV i];
            NADMV_idx=[NADMV_idx j];
        end


        load 'postIR_NTS_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_NTS = [paramsTemp_NTS; postIRstore(i).params(j,:)];
            HRstore_NTS = [HRstore_NTS postIRstore(i).HR(j)];
            COstore_NTS = [COstore_NTS postIRstore(i).CO(j)];
            crit_NTS = [crit_NTS postIRstore(i).crit(j,6)];
            A2_NTS = [A2_NTS postIRstore(i).Avals(j,2)];
            NA_NTS = [NA_NTS postIRstore(i).NA(j)];
            DMV_NTS = [DMV_NTS postIRstore(i).DMV(j)];
            symp_NTS = [symp_NTS postIRstore(i).symp(j)];
            preIRmdl_NTS = [preIRmdl_NTS i];
            NTS_idx=[NTS_idx j];
        end


        load 'postIR_ICN_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_ICN = [paramsTemp_ICN; postIRstore(i).params(j,:)];
            HRstore_ICN = [HRstore_ICN postIRstore(i).HR(j)];
            COstore_ICN = [COstore_ICN postIRstore(i).CO(j)];
            crit_ICN = [crit_ICN postIRstore(i).crit(j,6)];
            A2_ICN = [A2_ICN postIRstore(i).Avals(j,2)];
            NA_ICN = [NA_ICN postIRstore(i).NA(j)];
            DMV_ICN = [DMV_ICN postIRstore(i).DMV(j)];
            symp_ICN = [symp_ICN postIRstore(i).symp(j)];
            preIRmdl_ICN = [preIRmdl_ICN i];
            ICN_idx=[ICN_idx j];
        end


        load 'postIR_all_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_all = [paramsTemp_all; postIRstore(i).params(j,:)];
            HRstore_all = [HRstore_all postIRstore(i).HR(j)];
            COstore_all = [COstore_all postIRstore(i).CO(j)];
            crit_all = [crit_all postIRstore(i).crit(j,6)];
            A2_all = [A2_all postIRstore(i).Avals(j,2)];
            NA_all = [NA_all postIRstore(i).NA(j)];
            DMV_all = [DMV_all postIRstore(i).DMV(j)];
            symp_all = [symp_all postIRstore(i).symp(j)];
            preIRmdl_all = [preIRmdl_all i];
            all_idx=[all_idx j];
        end
    end
end

disp('HR store size')
size(HRstore_NTS)
size(HRstore_BR)
size(HRstore_cardiac)
size(HRstore_ICN)
size(HRstore_NADMV)
size(HRstore_all)
rng default


%% preIR regression
% sigIdx = [1 2 3 4 5 7 8 9 11 16 20 23 27 28 30 31 32 36 39 40 42 45 46 47 48 49 50 51];
% parameters = paramsTemp_preIR(idxStore,sigIdx);
% HRbase_preIR = HR_preIR(idxStore,zeroIdx);
% 
% X_norm = zscore(parameters);
% Mdl_preIR = fitrsvm(X_norm,HRbase_preIR,'KernelFunction','polynomial');
% 
% % Check regression fit
% predicted_preIR = predict(Mdl_preIR,X_norm);
% 
% figure;
% plot(predicted_preIR,HRbase_preIR,'bo')
% hold on
% plot([65 73],[65 73],'r-')
% title('Pre-I/R')
% xlabel('Predicted heart rate (bpm)')
% ylabel('Actual heart rate (bpm)')
% set(gca,'FontSize',16)
% saveas(gcf,'validation_preIR.png')
% 
% mse_final = mse(predicted_preIR,HRbase_preIR)
% shap_preIR = zeros(26,length(X_norm));
% % Calculate preIR Shapley values
% for i = 1:length(X_norm)
%     queryPoint = X_norm(i,:);
%     explainer_preIR = shapley(Mdl_preIR,X_norm,'QueryPoint',queryPoint,'UseParallel',true);
% 
%     shap_preIR(:,i) = table2array(explainer_preIR.ShapleyValues(:,2));
% end
% 
% shap_preIR_arr = shap_preIR;%table2array(shap_preIR);
% 
% % Plot shapley values
% figure;
% boxplot(shap_preIR_arr')
% title('Pre-I/R')
% ylabel('Feature value contribution')
% % ylim([-0.5 0.5])
% set(gca,'FontSize',16)
% set(gcf,'Position',[10 10 700 400])
% 
% % ratio of central vs peripheral shapley values for each accepted parameter
% % set
% temp_peri = zeros(1,size(shap_preIR_arr,1));
% temp_cent = zeros(1,size(shap_preIR_arr,1));
% ratio_peri_cent = zeros(size(shap_preIR_arr,1),1);
% for i = 1:size(shap_preIR_arr,2)
%     temp_peri(i) = abs(mean(shap_preIR_arr([1:13,26],i)));
%     temp_cent(i) = abs(mean(shap_preIR_arr([14:25],i)));
%     ratio_peri_cent(i,1) = temp_peri(i)./temp_cent(i);
% end
% 
% 
% figure;
% boxplot(ratio_peri_cent)
% ylabel('Ratio of peripheral/central contribution to heart rate');
% title('Pre-I/R');
% set(gca, 'FontSize', 16);
% 
% % Average Shapley value for central vs peripheral
% central_shap = mean(shap_preIR_arr(:,18:42),"all")
% peripheral_shap = mean(shap_preIR_arr(:,[1:17,51]),"all")
% 
% cols1 = [1:17, 51];
% cols2 = 18:42;
% 
% % Reshape the extracted columns into single columns
% data1 = reshape(shap_preIR_arr(:, cols1), [], 1);
% data2 = reshape(shap_preIR_arr(:, cols2), [], 1);
% 
% % Determine the length of the longer data set
% maxLength = max(length(data1), length(data2));
% 
% % Pad the shorter data set with NaNs
% if length(data1) < maxLength
%     data1 = [data1; NaN(maxLength - length(data1), 1)];
% elseif length(data2) < maxLength
%     data2 = [data2; NaN(maxLength - length(data2), 1)];
% end
% 
% % Combine the data into a new array
% data_combined = [data1, data2];
% 
% % Create a boxplot
% figure;
% boxplot(data_combined, 'Labels', {'Peripheral', 'Central'});
% ylabel('Contribution to heart rate');
% title('Pre-I/R');
% set(gca, 'FontSize', 16);
%% postIR regression
sigIdx = [1 2 3 4 5 7 8 9 11 16 20 23 27 28 30 31 32 36 39 40 42 45 46 47 48 49 50 51];
parameters = [paramsTemp_NTS; paramsTemp_ICN; paramsTemp_all; paramsTemp_BR; paramsTemp_cardiac; paramsTemp_NADMV];
% parameters = [paramsTemp_NTS(:,sigIdx); paramsTemp_ICN(:,sigIdx); paramsTemp_all(:,sigIdx); paramsTemp_BR(:,sigIdx); paramsTemp_cardiac(:,sigIdx); paramsTemp_NADMV(:,sigIdx)];
% HRbase_postIR = [HRstore_NTS'; HRstore_ICN'; HRstore_all'; HRstore_BR'; HRstore_cardiac'; HRstore_NADMV'];
% % crit = [crit_NTS; crit_ICN; crit_all; crit_BR; crit_cardiac; crit_NADMV];
% 
% X_norm = zscore(parameters);
% disp('parameters size')
% size(parameters)
% 
% % randomly sample X_norm to decrease computational expense
% rng default
% % sampleIdx = randi(length(X_norm),[499,1]);
% % X_norm_sample = X_norm(sampleIdx,:);
% % HRbase_postIR_sample = HRbase_postIR(sampleIdx,:);
% size(X_norm)
% size(HRbase_postIR)
% Mdl_postIR = fitrsvm(X_norm,HRbase_postIR,'KernelFunction','polynomial');
% 
% % Check regression fit
% predicted_postIR = predict(Mdl_postIR,X_norm);
% 
% figure;
% plot(predicted_postIR,HRbase_postIR,'bo')
% hold on
% plot([65 73],[65 73],'r-')
% title('Post-I/R')
% xlabel('Predicted heart rate (bpm)')
% ylabel('Actual heart rate (bpm)')
% set(gca,'FontSize',16)
% saveas(gcf,'validation_postIR_063024.png')
% 
% mse_final = mse(predicted_postIR,HRbase_postIR)
% 
% shap_postIR = zeros(26,length(X_norm));
% % Calculate postIR Shapley values
% for i = 1:length(X_norm)
%     queryPoint = X_norm(i,:);
%     explainer_postIR = shapley(Mdl_postIR,X_norm,'QueryPoint',queryPoint,'UseParallel',true);
% 
%     shap_postIR(:,i) = table2array(explainer_postIR.ShapleyValues(:,2));
% end
% 
% shap_postIR_arr = shap_postIR;%table2array(shap_postIR);

%% ratio of central vs peripheral shapley values for each accepted parameter
% set
% for i = 1:size(shap_postIR_arr,2)
%     temp_peri(i) = abs(mean(shap_postIR_arr([1:13,26],i)));
%     temp_cent(i) = abs(mean(shap_postIR_arr([14:25],i)));
%     ratio_peri_cent_post(i,1) = temp_peri(i)./temp_cent(i);
% end

%% Outlier removal-- ratios become significantly different pre and post IR
% % Calculate the IQR
% Q1 = quantile(ratio_peri_cent_post, 0.25);
% Q3 = quantile(ratio_peri_cent_post, 0.75);
% IQR = Q3 - Q1;
%
% % Calculate the mean
% mean_v = mean(ratio_peri_cent_post);
%
% % Determine the bounds for outliers
% lower_bound = mean_v - 1.5 * IQR;
% upper_bound = mean_v + 1.5 * IQR;
%
% % Remove outliers
% ratio_peri_cent_post = ratio_peri_cent_post(ratio_peri_cent_post >= lower_bound & ratio_peri_cent_post <= upper_bound);
%
% % Calculate the IQR
% Q1 = quantile(ratio_peri_cent, 0.25);
% Q3 = quantile(ratio_peri_cent, 0.75);
% IQR = Q3 - Q1;
%
% % Calculate the mean
% mean_v = mean(ratio_peri_cent);
%
% % Determine the bounds for outliers
% lower_bound = mean_v - 1.5 * IQR;
% upper_bound = mean_v + 1.5 * IQR;
%
% % Remove outliers
% ratio_peri_cent = ratio_peri_cent(ratio_peri_cent >= lower_bound & ratio_peri_cent <= upper_bound);
%%

% Determine the length of the longer data set
% maxLength = max(length(ratio_peri_cent_post), length(ratio_peri_cent));
% 
% % Pad the shorter data set with NaNs
% if length(ratio_peri_cent_post) < maxLength
%     ratio_peri_cent_post = [ratio_peri_cent_post; NaN(maxLength - length(ratio_peri_cent_post), 1)];
% elseif length(ratio_peri_cent) < maxLength
%     ratio_peri_cent = [ratio_peri_cent; NaN(maxLength - length(ratio_peri_cent), 1)];
% end

% save('central_peripheral.mat','ratio_peri_cent',"ratio_peri_cent_post",'shap_postIR_arr','shap_preIR_arr')

%% Plot
% boxplot
% load 'central_peripheral.mat'
% ratio_combined = [ratio_peri_cent, ratio_peri_cent_post];
% figure;
% boxplot(ratio_combined,'Labels',{'Pre-I/R','Post-I/R'})
% hold on
% plot([0 4],[1 1],'k--')
% ylabel('Peripheral/central contribution to heart rate');
% ylim([0 10])
% set(gca, 'FontSize', 16);
% set(gcf, 'Position',[10 10 400 500])
% saveas(gcf,'ratio_central_peripheral_063024.png')
% 
% % swarmchart
% x = [ones(1,length(ratio_peri_cent)); 2*ones(1,length(ratio_peri_cent_post))]';
% figure;
% swarmchart(x,ratio_combined)
% hold on
% plot([0.5 2.5],[1 1],'k--')
% ylabel('Peripheral/central contribution to heart rate');
% % ylim([0 2.5])
% set(gca, 'FontSize', 16);
% set(gcf, 'Position',[10 10 400 500])
% saveas(gcf,'ratio_central_peripheral_swarm_063024.png')
% 
% % t-test
% x = ratio_combined(:,1);
% y = ratio_combined(:,2);
% h = ttest(x,y,"Tail","both","Alpha",0.05)
% disp('if h=0, no evidence of observed effect')
% 
% %%
% % Assuming ratio_combined is a matrix where each row is a pair (before, after)
% ratios = ratio_combined(:, 1) ./ ratio_combined(:, 2);
% 
% % Perform a one-sample t-test to see if the mean of ratios is different from 1
% [h, p, ci, stats] = ttest(ratios, 1, 'Tail', 'both', 'Alpha', 0.05);
% 
% disp('if h=0, no evidence of observed effect');
% disp(['h = ', num2str(h)]);
% disp(['p-value = ', num2str(p)]);
% disp(['Confidence interval = ', num2str(ci')]);
% disp(['t-statistic = ', num2str(stats.tstat)]);

%% Plot BRS vs Shap value
% addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
A2_combined_raw = [A2_NTS, A2_ICN, A2_all, A2_BR, A2_cardiac, A2_NADMV];
A2_combined = -abs(A2_combined_raw);

NAactivity_combined = [NA_NTS, NA_ICN, NA_all, NA_BR, NA_cardiac, NA_NADMV];
DMVactivity_combined = [DMV_NTS, DMV_ICN, DMV_all, DMV_BR, DMV_cardiac, DMV_NADMV];
COstore_combined = [COstore_NTS, COstore_ICN, COstore_all, COstore_BR, COstore_cardiac, COstore_NADMV];
HRstore_combined = [HRstore_NTS, HRstore_ICN, HRstore_all, HRstore_BR, HRstore_cardiac, HRstore_NADMV];
symp_combined = [symp_NTS, symp_ICN, symp_all, symp_BR, symp_cardiac, symp_NADMV];
preIRmdl_combined = [preIRmdl_NTS, preIRmdl_ICN, preIRmdl_all, preIRmdl_BR, preIRmdl_cardiac, preIRmdl_NADMV];
idx_combined =  [NTS_idx, ICN_idx, all_idx, BR_idx, cardiac_idx, NADMV_idx];
% save('comb_postIR_accepted.mat','NAactivity_combined','DMVactivity_combined',"symp_combined",'A2_combined','COstore_combined','HRstore_combined','parameters','preIRmdl_combined','idx_combined')
load 'central_peripheral.mat'
load 'comb_postIR_accepted.mat'
load 'FilteredDistribution_061824.mat'
load 'FilteredDistributionIdx_061824.mat'

figure;
plot(abs(A2_combined),ratio_peri_cent_post,'ro')
hold on
plot(abs(A2Store(idxStore)),ratio_peri_cent(1:59),'bo')
xlabel('Baroreflex slope')
ylabel('Avg Shapley value')
set(gca,'FontSize',14)
saveas(gcf,'BRSvsShap.fig')
hold off;

figure;
hist(A2_combined)
saveas(gcf,'A2_hist.fig')

% BRS correlated to NA activity, DMV activity, symp activity, baroreceptor
% activity
figure;
plot(A2_combined,NAactivity_combined,'o')
xlabel('Baroreflex slope')
ylabel('NA activity')
set(gca,'FontSize',14)
saveas(gcf,'BRSvsNA.fig')

figure;
plot(A2_combined,DMVactivity_combined,'o')
xlabel('Baroreflex slope')
ylabel('DMV activity')
set(gca,'FontSize',14)
saveas(gcf,'BRSvsDMV.fig')

figure;
plot(A2_combined,symp_combined,'o')
xlabel('Baroreflex slope')
ylabel('Symp')
set(gca,'FontSize',14)
saveas(gcf,'BRSvsSymp.fig')

%% tsne
rng default
load 'neural_act_070424.mat'
Xdata = [DMVactivity_combined', NAactivity_combined', LSactivity,CPactivity,BRactivity];% , parameters(:,end),NTSactivity,  COstore_combined', ratio_peri_cent_post,symp_combined', sympActivity,
Xdata_clean = Xdata;
% Find rows with NaN values
nanRows = any(isnan(Xdata), 2);

% Get the indices of those rows
nanIndices = find(nanRows);

A2_clean = A2_combined;
A2_clean(:,nanIndices) = [];
Xdata_clean(nanIndices,:) = [];

% for i = 5:5:50
 i = 15;
% perplexity = 15 seems to be best looking, but all look pretty similar
y = tsne(Xdata,"Algorithm",'barneshut','Standardize',true,'Perplexity',i);

% Check t-SNE output dimensions
if size(y, 2) ~= 2
    error('t-SNE output dimensions are incorrect');
end

% Set colorData from A2_combined, ensuring dimensions match
colorData = abs(A2_clean(:)); % Ensure colorData is a column vector

% Check if colorData matches the number of samples
if size(colorData, 1) ~= size(Xdata_clean, 1)
    error('Number of samples in colorData does not match Xdata');
end

% Plotting
figure;
a = scatter(y(:, 1), y(:, 2), 15, colorData, 'filled'); % 15 is the marker size
colormap jet;  % Choose a colormap
colorbar;  

% Set axis labels and title for better visualization
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
set(gca,'FontSize',16)
% title(num2str(i))
saveas(gcf, 'tsne_coloredBRS_noBR_070524.png')
save('tsne_coloredBRS_noBR_070524.mat','y')
% end


