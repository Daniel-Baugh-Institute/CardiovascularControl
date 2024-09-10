% PreIR sigmoidal parameters compared for rejected and accepted
clear; 
close all;
% addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
% addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
% addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
% addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')
cd ..
my_dir = pwd;
addpath(genpath(my_dir))

% Parameters that change
date = '060824';
%% Plot parameter distributions for pre-IR and post-IR
load 'FilteredDistribution_061824.mat'
% A1Store, A2Store, A3Store, A4Store, DP, ECSP, HR, MAP , sampleSet (same as above), SP
matFilename = 'FilteredDistributionIdx_061824.mat';
load(matFilename) % gives 'idxStore','sampleSet'

BPvalues = -70:10:90; % neck chamber pressure
zeroIdx = find(BPvalues == 0);



% remove simulations with nan
% sampleSet = sampleSetAll;
% A1Store = A1all;
% A2Store = A2all;
% A3Store = A3all;
% A4Store = A4all;
% HR = HRall;
% MAP = MAPall;
allData = [A1Store A2Store A3Store A4Store HRstore(:,zeroIdx) MAPstore(:,zeroIdx)];


rejected_parameters = [A1Store, A2Store, A3Store, A4Store, HRstore(:,zeroIdx),MAPstore(:,zeroIdx)]; % 1:17, 42:46 cardiac parameters [1:4,4,5:8,9,9:13,14:17,17]
accepted_parameters = [A1Store(idxStore), A2Store(idxStore), A3Store(idxStore), A4Store(idxStore), HRstore(idxStore,zeroIdx),MAPstore(idxStore,zeroIdx)];
parameterNames =  {'A1 (s)','A2 (s/mm Hg)','A3 (mm Hg)','A4 (s)','Heart rate (bpm)','Mean arterial pressure (mm Hg)'};

saveFilename = ['preIR_paramDist_' date '.png'];

%% Load experimental data
patientData = {'P1','P2','P3','P4','P5'};

A1patient = zeros(length(patientData),1);
A2patient = zeros(length(patientData),1);
A3patient = zeros(length(patientData),1);
A4patient = zeros(length(patientData),1);
MAPbase = zeros(1,5);
HRbase = zeros(1,5);

% to get MAP curve, need to change initial guess in sigmoidRegression.m
for i = 1:length(patientData)
    load(patientData{i})
    data = eval(patientData{i});
    x_raw = data.raw_ECSP;
    y_raw = 60./data.raw_HR;
    [A1temp,A2temp,A3temp,A4temp] = sigmoidRegression(x_raw,y_raw);
    A1patient(i) = A1temp;
    A2patient(i) = A2temp;
    A3patient(i) = A3temp;
    A4patient(i) = A4temp;
    % Calculate baseline MAP and HR
    filename = [patientData{i} '_MAP.mat'];
    load(filename)
    varname = [patientData{i} '_MAP'];
    data = eval(varname);
    [A1_MAP,A2_MAP,A3_MAP,A4_MAP] = sigmoidRegression(data(:,1),data(:,2)); % MAP data regression
    MAPbase(i) = A3_MAP;
    HRbase(i) = (A1temp/(1+exp(A2temp*(MAPbase(i)-A3temp))))+A4temp;

    MAP_plot = (A1_MAP/(1+exp(A2_MAP*(data(:,1)-A3_MAP))))+A4_MAP;
    % plot(data(:,1),MAP_plot,'k--')
    % hold on
end

patientAvals = [A1patient, A2patient, A3patient, A4patient];
xvals = ones(length(patientData),1);

% Calculate 95% confidence interval for MAP and HR at baseline
HRbase = 60./HRbase;
MAPbase = MAPbase([1:3,5]); % remove patient 4 because of noise in data
HRstd = std(HRbase);
MAPstd = std(MAPbase);
HRupper = mean(HRbase) + 1.96* HRstd; %2.776*HRstd/sqrt(length(HRbase)); % assume t-distribution
HRlower = mean(HRbase) - 1.96* HRstd; %2.776*HRstd/sqrt(length(HRbase));
HR95 = [HRlower HRupper];
MAPupper = mean(MAPbase) + 1.96*MAPstd; %3.182*MAPstd/sqrt(length(MAPbase));
MAPlower = mean(MAPbase) - 1.96*MAPstd; %3.182*MAPstd/sqrt(length(MAPbase));
MAP95 = [MAPlower MAPupper];

% 95% confidence intervals 
A1range = [0.0316    0.4639];
A2range = [-0.3006    0];
A3range = [60.0981  128.6592];
A4range = [0.5480    1.0360];
HR95_Mastitskaya = [55 76];
MAP95_Mastitskaya = [69 130];
Arange = [A1range; A2range; A3range; A4range; HR95_Mastitskaya; MAP95_Mastitskaya];

%% plot
numParameters = 6;
grey = [0.8745    0.9098    0.8549];


figure;
tlo = tiledlayout(1,numParameters);
for i = 1:numParameters
    ax(i) = nexttile(tlo);
    x = [accepted_parameters(:, i); rejected_parameters(:, i)];
    g = [zeros(length(accepted_parameters(:, i)),1); ones(length(rejected_parameters(:, i)),1)];
    boxplot(ax(i),x,g,'Labels', {'Accepted', 'Rejected'});
    n = findobj(gcf,'tag','Outliers');
    for j = 1:numel(n)
        n(j).MarkerEdgeColor = grey;
    end
    if i  == 1
        hold on
        plot(xvals,patientAvals(:,i),'ro','MarkerFaceColor','r')
        % n = findobj(gcf,'tag','Outliers');
        % for j = 1:numel(n)
        %     n(j).MarkerEdgeColor = grey;
        % end
        legend({'Individual'; 'data'})
    elseif i < 5
        hold on
        plot(xvals,patientAvals(:,i),'ro','MarkerFaceColor','r')
        % n = findobj(gcf,'tag','Outliers');
        % for j = 1:numel(n)
        %     n(j).MarkerEdgeColor = grey;
        % end
    end
    ylabel(parameterNames{i});
    set(gca,'FontSize',24)
end

set(gcf,'Position',[100, 100, 2200,800])
saveas(gcf,saveFilename)

%% version showing all points
% SHOW CONVERTED MASTITSKAYA DATA POINTS FOR HR AND MAP INSTEAD OF
% SEREDYNSKI
rng default
max_jitter = 0.1;
min_jitter = 0;
x_jitter = (max_jitter - min_jitter).*randn(length(accepted_parameters),1) + min_jitter;
x_jitter2 = (max_jitter - min_jitter).*randn(length(rejected_parameters),1) + min_jitter;

figure;
tlo = tiledlayout(1,numParameters);
for i = 1:numParameters
    ax(i) = nexttile(tlo);
    x_accepted = ones(length(accepted_parameters(:,i)),1);
    x_rejected = ones(length(rejected_parameters(:,i)),1);
    if i == 2
        plot(x_rejected+x_jitter2,-abs(rejected_parameters(:, i)),'o','MarkerFaceColor',grey,'MarkerEdgeColor','k');
        hold on
        plot(x_accepted+x_jitter,-abs(accepted_parameters(:, i)),'bo','MarkerFaceColor','b','MarkerEdgeColor','k')
        plot([1.45 1.45],[A2range(1) A2range(2)],'k-','LineWidth',2.5)
    else

        plot(x_rejected+x_jitter2,abs(rejected_parameters(:, i)),'o','MarkerFaceColor',grey,'MarkerEdgeColor','k');
        hold on
        plot(x_accepted+x_jitter,abs(accepted_parameters(:, i)),'bo','MarkerFaceColor','b','MarkerEdgeColor','k')
        plot([1.45 1.45],[Arange(i,1) Arange(i,2)],'k-','LineWidth',2.5)
    end
    if i  == 1
        hold on
        plot(xvals,patientAvals(:,i),'r^','MarkerFaceColor','r','MarkerSize',14,'MarkerEdgeColor','k')
        % legend('Rejected','Accepted','Individual')
    elseif i < 5
        hold on
        plot(xvals,patientAvals(:,i),'r^','MarkerFaceColor','r','MarkerSize',14,'MarkerEdgeColor','k')
    elseif i == 5
        plot(xvals,HRbase,'r^','MarkerFaceColor','r','MarkerSize',14,'MarkerEdgeColor','k')
    elseif i == 6
        plot(xvals(1:4),MAPbase,'r^','MarkerFaceColor','r','MarkerSize',14,'MarkerEdgeColor','k')
    end

    % Set plotting range to 1.5*IQR
    iqr15 = 1.5*iqr(allData(:,i));
    per25 = quantile(allData(:,i),0.25);
    per75 = quantile(allData(:,i),0.75);
    lowerRange = per25 - iqr15;
    higherRange = per75 + iqr15;
    if i == 1 || i == 4
        lowerRange = 0;
    elseif i == 2
        higherRange = 0;
        lowerRange = -0.3;
    end

    ylim([lowerRange higherRange])

    ylabel(parameterNames{i});
    set(gca,'FontSize',28)
    set(gca,'XTick',[])
end

set(gcf,'Position',[100, 100, 2500,850])
saveFilename = ['preIR_params_' date '.png'];
saveas(gcf,saveFilename)

%% tSNE of accepted/rejected parameter sets
% Input variables: A1-A4, baseline HR, baseline MAP
% Color by experimental data, accepted, rejected
% rng default
% perplexity = 25;
% all_idx = 1:1:5300;
% rejected_idx = setdiff(all_idx,idxStore);
% X_individual = [A1patient([1:3,5]),A2patient([1:3,5]),A3patient([1:3,5]),A4patient([1:3,5]),HRbase([1:3,5])',MAPbase'];
% X_accepted = [A1Store(idxStore),A2Store(idxStore),A3Store(idxStore),A4Store(idxStore),HRstore(idxStore,zeroIdx),MAPstore(idxStore,zeroIdx)]; % 
% X_rejected = [A1Store(rejected_idx),A2Store(rejected_idx),A3Store(rejected_idx),A4Store(rejected_idx),HRstore(rejected_idx,zeroIdx),MAPstore(rejected_idx,zeroIdx)];
% % remove rows with nan
% rowsWithNaN = any(isnan(X_rejected), 2);
% X_rejected_cleaned = X_rejected(~rowsWithNaN, :);
% 
% X = [X_rejected_cleaned; X_accepted; X_individual];
% Y = tsne(X,'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1);
% 
% colors = [grey; 0 0 1; 1 0 0];
% g0 = repmat({'Individual'},1,4);
% g1 = repmat({'Accepted'},1,length(idxStore));
% g2 = repmat({'Rejected'},1,size(X_rejected_cleaned,1));
% labels =  cat(2,g2,g1,g0);
% sym = 'oo^';
% sz = [2 5 6];
% 
% 
% figure;
% a = gscatter(Y(:,1),Y(:,2),labels',colors,sym,sz);
% set(a(1),'MarkerFaceColor',colors(1,:))
% set(a(2),'MarkerFaceColor',colors(2,:))
% set(a(3),'MarkerFaceColor',colors(3,:))
% xlabel('t-SNE 1')
% ylabel('t-SNE 2')
% filename = ['preIR_sigmoid_tSNE_' date '.png']; %
% set(gca,'FontSize',16)
% saveas(gcf,filename)