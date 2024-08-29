%% Compare pre and post cardiac IR model parameters using PCA
% Purpose to show what is driving shift from pre to post IR physiology

clear; close all
%% Set up parallel pool
% myCluster = parcluster('local');
% myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
% myCluster.JobStorageLocation = getenv('TMPDIR');
% myPool = parpool(myCluster, myCluster.NumWorkers);
%%
% addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
% addpath 'C:\Users\mmgee\Downloads'
% addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
% addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')
my_dir = pwd;
addpath(genpath(my_dir))

%% Color by some BRS cutoff for those that met baseline postIR MAP, HR, baroreflex curve
% criteria. See if there's anything in preIR models that predicts that
% postIR BRS

%% Filter for postIR mdoels that meet quantitative data
load 'FilteredDistributionIdx_061824.mat' % variables: sampleSet, idxStore
numPatients = length(idxStore);
% BPvalues = -70:10:90; % neck chamber pressure
% zeroIdx = find(BPvalues == 0);
% altMdlName = {'cardiac', 'baroreceptors', 'NTS', 'NADMV', 'ICN', 'all'};
% 
% % Initialize cell arrays to store indices
% matchingIndices = cell(numPatients, 1);
% data_all = [];
% indices = [];
% Xdata = [];
% 
% % Loop through each of the 59 elements in postIRstore
% for j = 1:length(altMdlName)
%     load(['postIR_' altMdlName{j} '_06232024.mat'])
%     data = [];
%     for i = 1:numPatients
%         % Get the crit matrix for the current element
%         crit = postIRstore(i).crit;
% 
%         % Find rows where columns 1, 2, and 6 are 1, and column 7 is 0
%         % These criteria are SP, HR, and Osculati baroreflex curve
%         rows = find(crit(:, 1) == 1 & crit(:, 2) == 1 & crit(:, 6) == 1);% & crit(:, 7) == 0);
% 
%         % Store the matching rows in the cell array
%         if ~isempty(rows)
%             for k = 1:length(rows)
%             row = rows(k);
%             % Store the index and the required fields in the data array
%             indices = [indices; i, row, j];
%             data = [data; postIRstore(i).SP(row), postIRstore(i).DP(row), postIRstore(i).MAP(row), ...
%                     postIRstore(i).HR(row), postIRstore(i).CO(row), postIRstore(i).SDRR(row), ...
%                     postIRstore(i).RMSSD(row), postIRstore(i).Avals(row, 2), i];
%             data_all = [data_all; postIRstore(i).SP(row), postIRstore(i).DP(row), postIRstore(i).MAP(row), ...
%                     postIRstore(i).HR(row), postIRstore(i).CO(row), postIRstore(i).SDRR(row), ...
%                     postIRstore(i).RMSSD(row), postIRstore(i).Avals(row, 2), i];
%             end
%         end
% 
%     end
%     % Save the results
%     save(['postIRacceptedIdx_' altMdlName{j} '.mat'],'indices','data')
%     Xdata = [Xdata; data];
% end
% 
% % Add preIR data to tSNE matrix Xdata
% load 'FilteredDistribution_061824.mat'
% mdl_idx = 1:1:numPatients;
% Xdata = [Xdata; SPstore(idxStore,zeroIdx), DPstore(idxStore,zeroIdx), MAPstore(idxStore,zeroIdx), ...
%     HRstore(idxStore,zeroIdx), COstore(idxStore,zeroIdx), SDRR(idxStore), ...
%     RMSSD(idxStore), A2Store(idxStore), mdl_idx'];
% Xdata(:,8) = -abs(Xdata(:,8));
% matfilename = 'preAcceptedPostIR_tSNE_06232024.mat';
% save(matfilename,'Xdata')

%% pair plots
% columnNames = {'Systolic Pressure (mm Hg)', 'Diastolic Pressure (mm Hg)', 'MAP (mm Hg)', 'HR (bpm)', 'CO (mL/s)', 'SDRR (ms)', 'RMSSD (ms)', 'BRS'};
% numCols = length(columnNames);
% 
% % Iterate over each combination of columns
% for i = 1:numCols
%     for j = i+1:numCols
%         % Create a new figure for each combination
%         figure;
%         scatter(Xdata(:, i), Xdata(:, j), 'filled');
% 
%         % Set axis labels
%         xlabel(columnNames{i}, 'FontSize', 14);
%         ylabel(columnNames{j}, 'FontSize', 14);
% 
%         % Set title (optional)
%         title(['Plot of ', columnNames{i}, ' vs. ', columnNames{j}], 'FontSize', 14);
% 
%         % Adjust font size of axes
%         set(gca, 'FontSize', 14);
% 
%         % Enable grid
%         grid on;
%     end
% end
%% tSNE 
% rng default
% for i = 15%5:5:50
% perplexity = i; 
% % How do I map preIR to postIR? Do the clusters align in any way? Can be
% % that multiple preIR clusters map to one postIR cluster
% % 20 for postIR only, clusters are not super distinct, but one cluster
% % tends to have higher BRS
% % preIR: 3 clusters for 5, 2 for 15-20
% % but not separated by BRS
% % 45 for mixed pre and post IR
% Y = tsne(Xdata([558+1-numPatients:end],1:6),'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1); 
% %  preIR: ; postIR: 1:558-numPatients 
% % perplexity = 25;
% % Y25 = tsne(Xdata,'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1);
% % perplexity = 35;
% % Y35 = tsne(Xdata,'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1);
% % perplexity = 45;
% % Y45 = tsne(Xdata,'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1);
% 
% % matfilename = 'preAcceptedPostIR_tSNE_06232024.mat';
% % save(matfilename,'Xdata')
% % 
% % load(matfilename)
% 
% % generate labels
% g1 = repmat({'Pre-I/R'},numPatients,1);
% g0 = repmat({'Post-I/R'},size(Y,1)-numPatients,1);
% 
% % g1 = repmat({'Cardiac'},numCardiacmdls,1);
% % g3 = repmat({'Baroreceptors'},numBRmdls,1);
% % g4 = repmat({'NTS'},numNTSmdls,1);
% % g5 = repmat({'NA/DMV'},numNADMVmdls,1);
% % g6 = repmat({'ICN'},numICNmdls,1);
% % g7 = repmat({'All'},numAllmdls,1);
% % g1(idxStore_cardiac) = {'Cardiac'};
% % g1(idxStore_BR) = {'Baroreceptors'};
% % g1(idxStore_NTS) = {'NTS'};
% % g1(idxStore_NADMV) = {'NA/DMV'};
% % g1(idxStore_ICN) = {'ICN'};
% labels = cat(1,g0,g1);%,g3,g4,g5,g6,g7,g2);
% 
% 
% 
% % plot
% % close all;
% % brown = [ 0.7686    0.6000    0.4863];
% % green = [0.4660 0.6740 0.1880];
% % blue = [0 0 1];
% % grey = [0.9020    0.9020    0.9020];
% % greyOutline = [0.6510    0.6510    0.6510];
% % black = [0 0 0];
% % white = [1 1 1];
% % greyblue =[0.8510    0.9647    0.9882];
% % edgeClr = [greyOutline;  black; black; black;black;  black;black;  black];
% clr = ['r' 'b'];%[greyblue; 0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; green; 1 0 1; 0 1 1;  brown; blue]; 
% sym = 'po';%poooooop'; 
% % ms = 10;
% siz = [10 10];%[6 ms ms ms ms ms ms ms+2]; 
% 
% %  % Remove rows containing NaN values
% % Xdata_clean = Xdata(~any(isnan(Xdata), 2), :);
% 
% % Find rows with NaN values
% nanRows = any(isnan(Xdata), 2);
% 
% % Get the indices of those rows
% nanIndices = find(nanRows);
% 
% % Remove rows with NaN values
% Xdata_clean = Xdata(~nanRows, :);
% 
% % Display the number of rows removed
% rows_removed = size(Xdata, 1) - size(Xdata_clean, 1);
% disp(['Number of rows removed: ', num2str(rows_removed)]);
% % Extract the 8th column (BRS) for coloring
% % colorData = -abs(Xdata_clean(:, 8)); % 1:558-rows_removed-numPatients
% colorData = abs(Xdata(length(Xdata)+1-length(g1):end, 8)); % preIR ; postIR 1:length(Xdata)-rows_removed-numPatients
% 
% figure;
% size(Y)
% size(colorData)
% a = scatter(Y(:, 1), Y(:, 2),15, colorData, 'filled');%,labels,clr,sym,siz);  %  , 15 is the marker size
% colormap(jet);  % Choose a colormap
% colorbar;  % Add a colorbar
% % title(num2str(i))
% % a = gscatter(Y(:,1),Y(:,2),labels,clr,sym,siz);
% % set(a(1),'MarkerFaceColor',clr(1,:))
% % set(a(2),'MarkerFaceColor',clr(2,:))
% % set(a(3),'MarkerFaceColor',clr(3,:))
% % set(a(4),'MarkerFaceColor',clr(4,:))
% % set(a(5),'MarkerFaceColor',clr(5,:))
% % set(a(6),'MarkerFaceColor',clr(6,:))
% % set(a(7),'MarkerFaceColor',clr(7,:))
% % set(a(8),'MarkerFaceColor',clr(8,:))
% xlabel('t-SNE 1','FontSize',32)
% ylabel('t-SNE 2','FontSize',32)
% % l= findobj(gcf,'tag','legend'); set(l,'location','northeastoutside');
% 
% set(gca,'FontSize',24)
% set(gcf,'Position',[10 10 700 500])%[10 10 1250 900])
% filename = ['PreAcceptedPostIR_tSNE_06232024_' num2str(i) '.png'];
% saveas(gcf,filename)
% end



%% Mapping for pre to post IR
% find indices for preIR models when Y(1) and Y(2) within bound for cluster
% Find which postIR models are from which cluster of preIR models
% Plot postIR models by color according to preIR cluster

% g0_indices_preIR = find(Y(:,1) <= -64.7);
% g1_indices_preIR = setdiff(1:1:numPatients,g0_indices_preIR);
% save('preIR_grp_idx.mat','g0_indices_preIR','g1_indices_preIR')
% 
% % plot postIR models by color according to preIR cluster
% labels = repmat({'Group 1'},558-numPatients-rows_removed,1);
% 
% for i = 1:558-numPatients-rows_removed
%     if ~isempty(intersect(Xdata(i,9),g1_indices_preIR))
%         labels(i) = {'Group 2'};
%     end
% end 
% 
% 
% % postIR tsne
% for i = 15%5:5:50
% perplexity = i;
% % choose 35 to show somewhat distinct clusters that seem more enriched in
% % one group
% Y = tsne(Xdata(1:558-numPatients,1:8),'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1); 
% %  preIR: ; postIR: 1:558-numPatients 
% 
% 
% figure;
% a = gscatter(Y(:,1),Y(:,2),labels,clr,sym,siz);
% title(num2str(i))
% set(gca,'FontSize',24)
% set(gcf,'Position',[10 10 700 500])%[10 10 1250 900])
% filename = ['PostIR_tSNE_06232024_' num2str(i) '.png'];
% saveas(gcf,filename)
% end

%% plot by BRS
load 'neural_act_070424.mat'
load 'central_peripheral.mat'
load 'comb_postIR_accepted.mat'
load tsne_coloredBRS_allvars_070524.mat
load 'comb_postIR_accepted.mat'
 numPatients = 59;

 figure;
 colorData = abs(A2_combined');
 scatter(y(:,1),y(:,2),10,colorData,'filled');
 colormap(jet)
 colorbar('southoutside');
 xlabel('t-SNE Dimension 1')
ylabel('t-SNE Dimension 2')
 set(gca,'FontSize',24)
 set(gcf,'Position',[10 10 500 700])%[10 10 1250 900])

%% plot by coloration of preIR model 
% High BRS patient #'s: 4, 30, 41, 46, 23 


labels = repmat({'0'},499,1);
for i = 1:numPatients
    idx_temp = find(preIRmdl_combined == i);
    patient_label = {num2str(i)};
    labels(idx_temp) = patient_label;
end

% labels(nanIndices) = [];

figure;
a = gscatter(y(:,1),y(:,2),labels);%,clr,sym,siz);

set(gca,'FontSize',24)
legend('Location','northeastoutside','FontSize',10)
set(gcf,'Position',[10 10 700 500])%[10 10 1250 900])
filename = ['PostIR_tSNEbypreIR_06232024_' num2str(i) '.png'];
saveas(gcf,filename)

% Pre-I/R BRS of individuals in high BRS post-I/R gorup
load 'FilteredDistribution_061824.mat'
load 'FilteredDistributionIdx_061824.mat'
acceptedPreIR_A2 = A2Store(idxStore);
highBRS_A2 = acceptedPreIR_A2([4, 30, 41, 46, 23])
meanA2_preIRaccepted = mean(acceptedPreIR_A2)

% Post IR baroreceptor responsivity in individuals in high BRS post IR
% group
BRS_postIRidx = [222 232 226 172 134];
BR_responsitivity_postIR_highBRS = parameters(BRS_postIRidx,end)
BRS_responsivity_avg_postIR = mean(parameters(:,end))

%% plot by coloration of mode of adaptation
% Initialize a cell array
labels = cell(499, 1);
rng default
numMdls = [75 93 72 88 85 86];
labels(1:sum(numMdls(1))) = {'NTS'};
labels(sum(numMdls(1))+ 1:sum(numMdls(1:2))) = {'ICN'};
labels(sum(numMdls(1:2))+ 1:sum(numMdls(1:3))) = {'All'};
labels(sum(numMdls(1:3))+ 1:sum(numMdls(1:4))) = {'BR'};
labels(sum(numMdls(1:4))+ 1:sum(numMdls(1:5))) = {'Cardiac'};
labels(sum(numMdls(1:5))+ 1:sum(numMdls(1:6))) = {'NA/DMV'};

% labels(nanIndices) = [];

load 'neural_act_070424.mat'
load 'central_peripheral.mat'
load 'comb_postIR_accepted.mat'
Xdata = [DMVactivity_combined', NAactivity_combined', LSactivity,CPactivity,BRactivity, parameters(:,end)];% ,NTSactivity,  COstore_combined', ratio_peri_cent_post,symp_combined', sympActivity,
Xdata_clean = Xdata;
% Find rows with NaN values
nanRows = any(isnan(Xdata), 2);

 i = 20;
% perplexity = 15 seems to be best looking, but all look pretty similar
% y = tsne(Xdata,"Algorithm",'barneshut','Standardize',true,'Perplexity',i);

% colors
purple = [0.4940 0.1840 0.5560];
goldenrod = [0.9290 0.6940 0.1250];
green = [0.4660 0.6740 0.1880];
cyan = [0    1.0000    1.0000];
magenta = [1     0     1];
brown = [196/255 153/255 124/255];

colors = [green; cyan; brown; goldenrod; purple; magenta];
% [green; cyan; brown; goldenrod; purple; magenta];
load tsne_coloredBRS_allvars_070524.mat
figure;
a = gscatter(y(:,1),y(:,2),labels);%,clr,sym,siz);
% Update marker size for each group
for k = 1:length(a)
    a(k).MarkerSize = 16;
    a(k).Color = colors(k,:);
end
xlim([-65 40])
xlabel('t-SNE Dimension 1')
ylabel('t-SNE Dimension 2')
legend('Location','northeastoutside')
set(gca,'FontSize',24)
set(gcf,'Position',[10 10 800 500])%[10 10 1250 900])
filename = ['PostIR_tSNEbymode_06232024_' num2str(i) '.png'];
saveas(gcf,filename)

%% identify indices of postIR groups
load 'comb_postIR_accepted.mat'
idx = struct;
% group 1: good BRS
idx(1).grp_idx = find(y(:,1) < - 50);
% grp 2: middle cluster
idx(2).grp_idx = find(y(:,1) > -50 & y(:,1) < -15 & y(:,2) > -11);
%grp 3 bottom cluster
idx(3).grp_idx = find(y(:,1) < 0 & y(:,2) < -16);
% group 4: right cluster
idx(4).grp_idx = find(y(:,1) < -3);

% average BRS and vagal activity for each cluster
numClusters = 4;
BRS_mean = zeros(numClusters,1);
vagal_mean = zeros(numClusters,1);

for i = 1:numClusters
    indices = idx(i).grp_idx;
    BRS_mean(i) = mean(A2_combined(indices));
    vagal_mean(i) = ( mean(DMVactivity_combined(indices)) + mean(NAactivity_combined(indices)) + mean(LSactivity(indices)) + mean(BRactivity(indices)) + mean(CPactivity(indices)) )/5;
end

% Perform pairwise t-tests for BRS_mean
disp('Pairwise t-tests for BRS_mean:');
for i = 1:numClusters
    indices_i = idx(i).grp_idx;
    for j = i+1:numClusters
        indices_j = idx(j).grp_idx;
        [h, p] = ttest2(A2_combined(indices_i), A2_combined(indices_j));
        fprintf('BRS_mean: Cluster %d vs Cluster %d: p-value = %.4f\n', i, j, p);
        if h == 0
            disp('No significant difference between means');
        else
            disp('Significant difference between means');
        end
    end
end

% Perform pairwise t-tests for vagal_mean
disp('Pairwise t-tests for vagal_mean:');
for i = 1:numClusters
    indices_i = idx(i).grp_idx;
    BRgain_i = ( DMVactivity_combined(indices_i) + NAactivity_combined(indices_i) + LSactivity(indices_i)' + BRactivity(indices_i)' + CPactivity(indices_i)' )./5;
    for j = i+1:numClusters
        indices_j = idx(j).grp_idx;
        BRgain_j = (DMVactivity_combined(indices_j) + NAactivity_combined(indices_j) + LSactivity(indices_j)' + BRactivity(indices_j)' + CPactivity(indices_j)') ./5;
        [h, p] = ttest2(BRgain_i, BRgain_j);
        fprintf('Vagal_mean: Cluster %d vs Cluster %d: p-value = %.4f\n', i, j, p);
        if h == 0
            disp('No significant difference between means');
        else
            disp('Significant difference between means');
        end
    end
end

% Average baroreceptor gain for each subgroup to see if they are significantly different 
BR_gain_mean = zeros(1,numClusters);
for i = 1:numClusters
    indices = idx(i).grp_idx;
    BR_gain_mean(i) = mean(parameters(indices,end));
end

disp('Pairwise t-tests for barorecptor gain:');
for i = 1:numClusters
    indices_i = idx(i).grp_idx;
    BRgain_i = parameters(indices_i,end);
    for j = i+1:numClusters
        indices_j = idx(j).grp_idx;
        BRgain_j = parameters(indices_j,end);
        [h, p] = ttest2(BRgain_i, BRgain_j);
        fprintf('Baroreceptor gain: Cluster %d vs Cluster %d: p-value = %.4f\n', i, j, p);
        if h == 0
            disp('No significant difference between means');
        else
            disp('Significant difference between means');
        end
    end
end

%% Plot comparison of baroreceptor responsivity curves for pre and postIR
% post IR
postIR_mean = mean(parameters(:,51));
postIR_std = std(parameters(:,51));
postIR_lb = postIR_mean - postIR_std;
postIR_ub = postIR_mean + postIR_std;

% Load preIR data
load 'FilteredDistributionIdx_061824.mat' % variables: sampleSet, idxStore
preIR_accepted = sampleSet(idxStore,51);
preIR_mean = mean(preIR_accepted);
preIR_std = std(preIR_accepted);
preIR_lb = preIR_mean - preIR_std;
preIR_ub = preIR_mean + preIR_std;

fmin = 2.52; % Hz
fmax = 47.78; % Hz
p_prime = 60:1:120;
pn = 92; % mmHg

% Calc preIR baroreceptor response curve
fbr_mean_pre = (fmin+fmax.*exp((p_prime-pn)./preIR_mean))./(1+exp((p_prime-pn)/preIR_mean));
fbr_lower_pre = (fmin+fmax.*exp((p_prime-pn)./preIR_lb))./(1+exp((p_prime-pn)/preIR_lb));
fbr_upper_pre = (fmin+fmax.*exp((p_prime-pn)./preIR_ub))./(1+exp((p_prime-pn)/preIR_ub));

% Calc postIR baroreceptor response curve
fbr_mean_post = (fmin+fmax.*exp((p_prime-pn)./postIR_mean))./(1+exp((p_prime-pn)/postIR_mean));
fbr_lower_post = (fmin+fmax.*exp((p_prime-pn)./postIR_lb))./(1+exp((p_prime-pn)/postIR_lb));
fbr_upper_post = (fmin+fmax.*exp((p_prime-pn)./postIR_ub))./(1+exp((p_prime-pn)/postIR_ub));


% Plot preIR baroreceptor response curve
figure;
hold on;

% Shaded area for preIR standard deviation
fill([p_prime, fliplr(p_prime)], [fbr_lower_pre, fliplr(fbr_upper_pre)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');

% Mean line for preIR
plot(p_prime, fbr_mean_pre, 'b-', 'LineWidth', 2, 'DisplayName', 'Pre-IR');

% Shaded area for postIR standard deviation
fill([p_prime, fliplr(p_prime)], [fbr_lower_post, fliplr(fbr_upper_post)],[0.9290 0.6940 0.1250], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');

% Mean line for postIR
plot(p_prime, fbr_mean_post, '-','color',[0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'Post-IR');

% Customize plot
xlabel('Arterial pressure (mm Hg)', 'FontSize', 14);
ylabel('Firing rate (Hz)', 'FontSize', 14);
legend('Location', 'Best');
set(gca, 'FontSize', 14);
grid on;
hold off;
%% plot baroreflex curves from group 4 (low vagal activity, med BRS) to show
% that vagal stim can improve BRS and that BR stimulation does not improve
% BRS
% mdlName = 'ICN_model_v15_Mastitskaya2012_control_r2020b';
% BPvalues = 0;%[-70.0	-50.0	-25.0	0.0	8.0	13.2	18.5	23.7	30.0	37.5	50.0	75.0	90.0]; % neck chamber pressure
% 
% 
% % choose random parameter set from group
% rng(73) % seed 12345
% rand_mdl = randi(91,1); %215 for grp 4
% mdl_parameters = parameters(idx(2).grp_idx(1:91),:);
% %435 is an NADMV model
% 
% for i = 1:91
%     filename = ['highVagal_lowBRS_grp_VNS_10Hz_' num2str(i)];
% %ICN
% ICNparams = mdl_parameters(i,1:17);
% 
% % NA
% NAparams = mdl_parameters(i,30:41);
% 
% % NTS
% NTSparams = mdl_parameters(i,18:29);
% 
% % heart
% heart_params = mdl_parameters(i,[42,45:50]);
% 
% % BR
% ka = mdl_parameters(i,51);
% 
% baseline = 1; % calculate hemodynamic and neural metrics
% 
% 
% 
% [mat_filename, HRvalidation, MAP_EI_mdl, SP_store, DP_store, CO_store, ...
%     NA_activity_postIR, DMV_activity_postIR, symp_activity_postIR, ...
%     SDRR_store, RMSSD_store, pNN50_store] = baroreflex_curve(mdlName, ...
%     BPvalues,filename,ICNparams,NAparams,NTSparams,heart_params,ka,baseline);
% 
% HR_10Hz(i) = HRvalidation;
% end
% save('HR_highVagal_lowBRS_grp_VNS_10Hz.mat',"HR_10Hz")
%%
% load 'highVagal_lowBRS_VNS_10Hz_seed73.mat'
% load 'lowVagal_medBRS_VNS_4Hz.mat'
% ECSP = MAP_EI_mdl + BPvalues;
% RR = 1000*60./HR';
% figure
% [A1,A2,A3,A4] = sigmoidRegression(ECSP',RR)
% figure;
% plot(BPvalues,1000*60./HR)
% hold on;
% xlabel('Neck chamber pressure (mm Hg)')
% ylabel('RR interval (ms)')
% set(gca,'FontSize',14)

%% close_system
% delete(myPool);
% exit
% Plot baroreflex curves from group 2 (high vagal activity, low BRS to show
% that BR stimulation can improve BRS and that VNS does not improve BRS)