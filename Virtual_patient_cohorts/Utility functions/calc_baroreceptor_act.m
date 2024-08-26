%% Compare the pre- and post-IR relative influence of central (NTS, NA, DMV) and peripheral (ICN, baroreceptors) on heart rate
% Approach: Regression where the inputs are the parameters and the output
% is heart rate. Average the pre-IR central parameter Shapley values and
% compare to the pre-IR peripheral Shapley values. Then compare how this
% ratio changes pre- to post-IR
% clear;
close all;
%% Load data
addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')
%%
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
crit_BR = [];
A2_BR = [];
NA_BR = [];
DMV_BR = [];
symp_BR = [];

paramsTemp_ICN = [];
HRstore_ICN = [];
crit_ICN = [];
A2_ICN = [];
NA_ICN = [];
DMV_ICN = [];
symp_ICN = [];

paramsTemp_cardiac = [];
HRstore_cardiac = [];
crit_cardiac = [];
A2_cardiac = [];
NA_cardiac = [];
DMV_cardiac = [];
symp_cardiac = [];

paramsTemp_NADMV = [];
HRstore_NADMV = [];
crit_NADMV = [];
A2_NADMV = [];
NA_NADMV = [];
DMV_NADMV = [];
symp_NADMV = [];

paramsTemp_NTS = [];
HRstore_NTS = [];
crit_NTS = [];
A2_NTS = [];
NA_NTS = [];
DMV_NTS = [];
symp_NTS = [];

paramsTemp_all = [];
HRstore_all = [];
crit_all = [];
A2_all = [];
NA_all = [];
DMV_all = [];
symp_all = [];


for i = 1:numPatients
    for j = 1:numSamples
        load 'postIR_baroreceptors_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_BR = [paramsTemp_BR; postIRstore(i).params(j,:)];
        end


        load 'postIR_cardiac_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_cardiac = [paramsTemp_cardiac; postIRstore(i).params(j,:)];
        end


        load 'postIR_NADMV_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_NADMV = [paramsTemp_NADMV; postIRstore(i).params(j,:)];
        end


        load 'postIR_NTS_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_NTS = [paramsTemp_NTS; postIRstore(i).params(j,:)];
        end


        load 'postIR_ICN_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_ICN = [paramsTemp_ICN; postIRstore(i).params(j,:)];
        end


        load 'postIR_all_06232024.mat'
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            paramsTemp_all = [paramsTemp_all; postIRstore(i).params(j,:)];
        end
    end
end



%% postIR regression
% sigIdx = [1:5,7:12,16,20,23,27,28,30,31,32,36,38:40,42,42,51];
parameters = [paramsTemp_NTS; paramsTemp_ICN; paramsTemp_all; paramsTemp_BR; paramsTemp_cardiac; paramsTemp_NADMV];

% Run simulations in parallel
filename = 'neural_act_070424';
[NAactivity,DMVactivity,sympActivity,BRactivity,LSactivity,CPactivity,NTSactivity] = calcVagalActivity(parameters,filename);

% Save BR, CP, LS, NTS activity