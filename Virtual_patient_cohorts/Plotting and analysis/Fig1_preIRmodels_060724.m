%% Fig 1 plot
close all; clear;
% Add file paths
% addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
% addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
% addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')
my_dir = pwd;
addpath(genpath(my_dir))
%% plot unfiltered data
% load data
load 'all[25_25]500_06182024.mat' 

% plot formatting
lw = 2;
fs = 18;
gray = [0.6 0.6 0.6 0.6];

% Plot
RR_plot = 60./HR;
figure;
% for i = 1:length(A1Store)
i = 1;
    if i == 1
        plot(ECSP(1,:),RR_plot(1,:),'LineWidth',lw/3,'Color',gray)
        i = i + 1;
        hold on
    end
        plot(ECSP',RR_plot','LineWidth',lw/3,'Color',gray,'HandleVisibility','off')

% end

% ylim([0.6 1.2])
xlabel('Estimated carotid sinus pressure (mm Hg)')
ylabel('RR interval (s)') %
legend('Individual parameter set','Location','northeast') % 'No post-IR model','Post-IR model',
set(gca,'FontSize',fs)
hold off
saveas(gcf,'UnfilteredPreIRModels.png')


%% plot filtered data
% load preIR sampled models for 6000 samples
load 'FilteredDistribution_061824.mat'
HR = HRstore;
ECSP = ECSPstore;
CO = COstore;
DMVactivity = DMVactivityStore';
DP = DPstore;
SP = SPstore;
MAP = MAPstore;
NAactivity = NAactivityStore';
RMSSDstore = RMSSD';
SDRRstore = SDRR';
pNN50store = PNN50';
sampleSet = sampleSetStore;


% Identify sample sets where there was a nan value for HR/ECSP because A1-A4
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

matFilename = 'FilteredDistribution_060424.mat';
save(matFilename, 'A1Store', 'A2Store', 'A3Store', 'A4Store', 'COstore','DMVactivityStore', ...
    'DPstore', 'ECSPstore', 'HRstore','MAPstore','NAactivityStore','RMSSD','PNN50', ...
    'SDRR', 'SPstore','sampleSetStore')


%% find indices of sample sets that are within the range of observed
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
                                if RR_plot(i,1) > 0.6 && RR_plot(i,end) < 1.2
                                    idxStore = [idxStore i];
                                end
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
save('FilteredDistributionIdx_060424.mat','idxStore','sampleSet')
%% Plot accepted models
% load patient data
patientData = {'P1','P2','P3','P4','P5'};

for i = 1:length(patientData)
    load(patientData{i})
end

% Plot formatting
plot_rows = 3;
plot_cols = 1;
lw = 2;
ms = 6;
fs = 18;
gray = [0.9 0.9 0.9 0.3];

figure;
gray = [0.5 0.5 0.5 1];
color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880; 0.3020    0.7451    0.9333];
% subplot(plot_rows,plot_cols,2)
for i = 1:length(patientData)
    data = eval(patientData{i});

    % patient data
    x_raw = data.raw_ECSP;
    y_raw = 60./data.raw_HR;
    plot(x_raw,y_raw,'ko','MarkerSize',ms,'MarkerFaceColor',color(i,:),'HandleVisibility','off')
    hold on
    x_fit = data.fit_ECSP;
    y_fit = 60./data.fit_HR;
    plot(x_fit,y_fit,'--','LineWidth',lw,'Color',color(i,:));
    hold on
end

% load 'sprenkle.mat'
% plot(sprenkle(:,1),sprenkle(:,2),'ko','MarkerSize',ms,'MarkerFaceColor',color(i,:),'HandleVisibility','off')

ylim([0.6 1.2])%([0 4])%
xlabel('Estimated carotid sinus pressure (mm Hg)')
ylabel('RR interval (s)') %
legend('Individual 1','Individual 2','Individual 3','Individual 4','Individual 5','Location','northeast') % 'Individual 6',
% ttl = title('A');
% ttl.Units = 'Normalize';
% ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
% ttl.HorizontalAlignment = 'left';
set(gca,'FontSize',fs)
hold off
saveas(gcf,'SeredynskiData.png')

figure;
for i = 1:length(idxStore)
    if i == 1
        plot(ECSP(idxStore(i),:),RR_plot(idxStore(i),:),'LineWidth',lw/4,'Color',gray)
        hold on
    else
        plot(ECSP(idxStore(i),:),RR_plot(idxStore(i),:),'LineWidth',lw/4,'Color',gray,'HandleVisibility','off')
    end
end

% hasPostIRmdl = [6,9,10,14,15,18,20,24,28,29,32,36,41];
% for i = 1:length(hasPostIRmdl)
%     plot(ECSP(idxStore(hasPostIRmdl(i)),:),RR_plot(idxStore(hasPostIRmdl(i)),:),'LineWidth',lw/4,'Color','r')
% end

ylim([0.6 1.2])%([0 4])%
xlabel('Estimated carotid sinus pressure (mm Hg)')
ylabel('RR interval (s)') %
legend('Individual parameter set','Location','northeast') % 'No post-IR model','Post-IR model',
% ttl = title('B');
% ttl.Units = 'Normalize';
% ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
% ttl.HorizontalAlignment = 'left';
set(gca,'FontSize',fs)
hold off
saveas(gcf,'AcceptedPreIRModels.png')

% matFilename = ['FilteredDistribution032324_test.mat'];
% save(matFilename,'idxStore','sampleSet','MAP','HR')
disp('Number of accepted parameter sets')
disp(length(idxStore))

% Calculate NA and DMV activity for 29 accepted preIR models
% load 'FilteredDistribution030824.mat' %sampleSet, idxStore
%
% parameters = sampleSet(idxStore,:);
% filename = 'preIRvagalActivity';
% [NAactivity,DMVactivity,sympActivity] = calcVagalActivity(parameters,filename);


%% PCA
% Z = zscore(sampleSet);%(:,[1,21,41])); % columns centered to zero mean and SD is 1
% rejectedIdx = setdiff(1:1:700,idxStore);
% [coeff,score,latent] = pca(Z);
% figure;
% plot(score(idxStore,1),score(idxStore,2),'bo','MarkerFaceColor','b');%,'varlabels',{'v_1','v_2','v_3','v_4'});
% hold on
% plot(score(rejectedIdx,1),score(rejectedIdx,2),'o','MarkerFaceColor',[0.7 0.7 0.7]);

%% Regression
% numParamSets = size(sampleSet,1);
% X = zscore(sampleSet);
% respVar = zeros(numParamSets,1);%(223,numParamSets,1); %
% respVar(idxStore) = 1;%222;%
% y = respVar; % accepted = 1, rejected = 0;
% mdl1 = fitglm(X,y,'linear')
% ,'VarNames',{'PN_{NA} f_{min}','PN_{NA} f_{max}','PN_{NA} f_{mid}','PN_{NA} k',...
%     'LCN f_{min}','LCN f_{max}','LCN f_{mid}','LCN k',...
%     'PN_{DMV} f_{min}','PN_{DMV} f_{max}','PN_{DMV} f_{mid}','PN_{DMV} k',...
%     'PN_{NA} f_{min}','PN_{NA} f_{max}','PN_{NA} f_{mid}','PN_{NA} k',...
%     'NTS_{BR} f_{min}','NTS_{BR} f_{max}','NTS_{BR} f_{mid}','NTS_{BR} k',...
%     'NTS_{CP} f_{min}','NTS_{CP} f_{max}','NTS_{CP} f_{mid}','NTS_{CP} k',...
%     'NTS_{LS} f_{min}','NTS_{LS} f_{max}','NTS_{LS} f_{mid}','NTS_{LS} k',...
%     'NA f_{min}','NA f_{max}','NA f_{mid}','NA k',...
%     'NActr f_{min}','NActr f_{max}','NActr f_{mid}','NActr k',...
%     'DMV f_{min}','DMV f_{max}','DMV f_{mid}','DMV k',...
%     'C_{la}','P_{0lv}','k_{Elv}','V_{ulv}','fes_{inf}','k_a'})
% sigIdx = find(mdl1.Coefficients.pValue <= 0.05) - 1
%
% writetable(mdl1.Coefficients,'regression_allVars.csv')
%
% mdl2 = fitglm(X(:,sigIdx(2:end,1)),y,'linear')
% writetable(mdl2.Coefficients,'regression_sigVars.csv')


% Significant variables A1-A4 range: 1, 21, 41 PNNA fmin, NTS BR k, DMV k, intercept
% Significant variables: 2,3,19,20,23,32 PNNA fmax (2), PNNA fmid (3), NTS BR fmax 19, NTS BR fmid 20, NTS CP
% fmax 23, NA fmid (32)
%     if max(RR_plot(i,:)) < 1.17
%     if RR_plot(i,1) < 0.96
%     if RR_plot(i,end) > 0.85
%     if min(RR_plot(i,:)) > 0.65

% Significant variables 1,2, 3, 20, 30, 32
% Both of the above constraints together

% All of above plus Mastitskaya constraints on HR and MAP: 9, 26 (fmin PN
% DMV and fmin NTS LS subgroup

%% tSNE of parameters
% z score data


% numParamSets = rows;
% labels = repmat({'Rejected'},numParamSets,1);%(223,numParamSets,1); %
% labels(idxStore) = {'Accepted'};%222;%
% % umap_input = [Z(:,[18:34]) labels];
%
% rng default
% perplexity = 35;
% % end_col = size(umap_input,2);
% Y = tsne(sampleSet(:,sigIdx(2:end)),'Algorithm','exact','Distance','cityblock','Standardize',1);
% % [reduction, umap, clusterIdentifiers, extras]=run_umap(umap_input,'marker_size',10,'save_output','true','label_column',end_col)
%
% ms = 10;
% subplot(plot_rows,plot_cols,3)
% gscatter(Y(:,1),Y(:,2),labels,[0.7 0.7 0.7;0 0 1]);
% l= findobj(gcf,'tag','legend'); set(l,'location','Northwest');
% xlabel('t-SNE 1')
% ylabel('t-SNE 2')
% ttl = title('C');
% ttl.Units = 'Normalize';
% ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
% ttl.HorizontalAlignment = 'left';
% filename = 'Fig1_preIRmodels.png'; %
% set(gca,'FontSize',fs)
% set(gcf,'Position',[10,10 500,1200])
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

