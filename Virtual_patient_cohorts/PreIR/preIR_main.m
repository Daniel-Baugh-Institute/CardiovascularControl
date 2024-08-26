clear; close all;
%%
addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')

%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);

%% Fit Seredynski ensemble data
% load data
matFileNames = {'P1','P2','P3','P4','P5'};
for i = 1:length(matFileNames)
    fileName = [matFileNames{i}, '.mat'];
    load(fileName)
end

% concatenate data from patients to combine into an ensemble dataset for
% fitting
ECSP = [P1.raw_ECSP; P2.raw_ECSP; P3.raw_ECSP; P4.raw_ECSP; P5.raw_ECSP];
HR = [P1.raw_HR; P2.raw_HR; P3.raw_HR; P4.raw_HR; P5.raw_HR];
RR = 60./HR;
[A1,A2,A3,A4] = sigmoidRegression(ECSP,RR)

%% Fit individual data to determine range of acceptable baroreflex behavior
for i = 1:length(matFileNames)
    ECSP = eval([matFileNames{i} '.raw_ECSP']);
    RR = 60./eval([matFileNames{i} '.raw_HR']);
    [A1,A2,A3,A4] = sigmoidRegression(ECSP,RR);
    A1Store(i) = A1;
    A2Store(i) = A2;
    A3Store(i) = A3;
    A4Store(i) = A4;
end

A1mean = mean(A1Store);
A2mean = mean(A2Store);
A3mean = mean(A3Store);
A4mean = mean(A4Store);

A1std = std(A1Store);
A2std = std(A2Store);
A3std = std(A3Store);
A4std = std(A4Store);

A1CI95 = 1.96.* A1std;
A2CI95 = 1.96.* A2std;
A3CI95 = 1.96.* A3std;
A4CI95 = 1.96.* A4std;

A1range = [A1mean - A1CI95 A1mean + A1CI95]
A2range = [A2mean - A2CI95 A2mean + A2CI95] 
disp('This must be < 0 so that slope is > 0')
A2range = [A2mean - A2CI95 0]
A3range = [A3mean - A3CI95 A3mean + A3CI95]
A4range = [A4mean - A4CI95 A4mean + A4CI95]

%% Baroreflex curve from Gee 2023 model compared to Seredynski data
mdlName = 'ICN_model_v15_Mastitskaya2012_control';
BPvalues = -70:5:90;
t = char(datetime);
filename_prefix = 'baroreflex_Gee2023';
filename = filename_prefix;%[filename_prefix t];

% control PI parameters
PNDMVdelay = 0.3;
ICNparams = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 PNDMVdelay 3.329861 2.661685 5.642977 0.066794];

% NA/DMV
NAparams = [4.88, 15.78, 59.83, 23, 0.61, 11, 12.81, 7, 2.5901, 6.66, 42.91, 33.5]; %  BR, CPR, LSR

% NTS
NTSparams = [0.30, 21.50, 37.07, 21, 0.45, 28.33, 10.2, 7, 2.75, 31.57, 11.13, 2]; %  BR, CPR, LSR; fmin, fmax, fmid, k

% heart_params
Cla = 19.23;
P0lv = 1.5;
kElv = 0.014;
fes_inf = 2.1;
Vulv = 16.77;
heart_params = [Cla,P0lv,kElv,Vulv,fes_inf]; % PI
ka = 11.758;

mat_filename = baroreflex_curve(mdlName,BPvalues,filename,ICNparams,NAparams,NTSparams,heart_params,ka);
load(mat_filename)
MAP = baselineMAP;
% compare_baroreflex_Seredynski(MAP,BPvalues,HR)

%% Fit baroreflex curve
parameters = ICNparams(1:2);
numSampleSets = 1000; % minimum 5
variationFactor = 2;
fitParams = [A1,A2,A3,A4];
paramMatrix = paramEstimationBaroreflexCurve(parameters,mdlName,numSampleSets,variationFactor,fitParams);

%% Exit code
close_system
delete(myPool);
exit
%% Plots to find acceptable range of variation for ka (other parameters too?)