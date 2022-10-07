% Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
% Michelle Gee
% October 3, 2022

% Script to produce Figure 7

%% Rajendran ICN pig data
% Data set from Pennsieve: https://app.pennsieve.io/N:organization:618e8dd9-f8d2-4dc4-9abb-c6aaab2e78a0/datasets/N:dataset:0d9454ca-43d9-4fdb-ac79-01c3574e8565/files/N:collection:063f029d-ebf1-49fc-803c-1c618f314ba4
% Protocol: https://www.protocols.io/view/pig-icn-recording-2jugcnw?step=4
clear; close all;

% load data
load('Pig013_ICNS15_Matlab.mat') % Rajendran data

% RR interval calculation and plot
% time intervals for baseline, stim, recovery (see excel spreadsheet event log provided on Pennsieve)
lcvBaseTime = [3052.2 3112.2];
lcvStimTime = [3112.2 3171.1];
lcvRecovTime = [3171.1 3231.1];
rcvBaseTime = [3457.1 3517.1];
rcvStimTime = [3517.1 3576.0];
rcvRecovTime = [3576.0 3636.0];

% Extract ECG period of interest
res = ECG.interval;
times = res *ones(1,length(ECG.values));
tVec = cumsum(times) + ECG.start;

low = find(tVec >= lcvBaseTime(1));
high = find(tVec <= lcvRecovTime(2));
Idx = intersect(low,high);
baseECG = ECG.values(Idx);

% Calculate RR interval
[pks,locs] = findpeaks(baseECG,tVec(Idx),'MinPeakProminence',2);
RRint = diff(locs);
RRtimes = locs;
RRtimes(end) = [];
RRIdxBase = find(RRtimes >= lcvBaseTime(1) & RRtimes <= lcvBaseTime(2));
RRIdxStim = find(RRtimes >= lcvStimTime(1) & RRtimes <= lcvStimTime(2));
RRIdxRecov = find(RRtimes >= lcvRecovTime(1) & RRtimes <= lcvRecovTime(2));

% Convert RR interval to HR
HR = 60./RRint;

% Calculate RR change from baseline
lowBase = find(tVec >= lcvBaseTime(1));
highBase = find(tVec <= lcvBaseTime(2));
IdxBase = intersect(lowBase,highBase);
baseRR = mean(RRint(RRIdxBase));
stimRR = mean(RRint(RRIdxStim));
recovRR = mean(RRint(RRIdxRecov));

% Convert RR intervals to HR
baseHR = 60./baseRR;
stimHR = 60./stimRR;
recovHR = 60./recovRR;

%Percent change in RR interval from base
stimDel = (stimRR - baseRR) / baseRR
recovDel = (recovRR - baseRR) / baseRR

HRstim = (stimHR - baseHR) / baseHR
HRrecov = (recovHR - baseHR) / baseHR

%% simulation
% Define simulation input
Params = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 3.329861 2.661685 5.642977 0.066794];
kRSA = 0.5;

mdlName = 'ICN_model_v15_VNS';
simIn = Simulink.SimulationInput(mdlName);
simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(Params(2)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(Params(3)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(Params(4)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(Params(5)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(Params(6)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(Params(7)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(Params(8)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(Params(9)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(Params(10)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(Params(11)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(Params(12)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(13)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(14)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(Params(15)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)), ...
    [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(kRSA));

% run simulation
simOut = sim(simIn);

%% Calculate baseline, stim, recovery RR interval
simTime         = simOut.time;
tbase = [120 180]; % baseline
tstim = [180 240]; % stimulation
trecov = [240 300]; % recovery

% convert intervals to model indices
tbaseIdx        = find(simTime>= tbase(1)  & simTime <= tbase(2));
tstimIdx        = find(simTime>= tstim(1)  & simTime <= tstim(2));
trecovIdx        = find(simTime>= trecov(1)  & simTime <= trecov(2));

time = simOut.time(tbaseIdx);
phi = simOut.Phi(tbaseIdx);

% RR interval calculation for baseline interval
% Determine indices for when phi=0 (beginning and end of heart beat)
for i = length(phi):-1:2
    if phi(i) == phi(i-1)
        phi(i-1) = [];
        time(i-1) = [];
    end
end

RRidx = find(~phi);
RRtimes = time(RRidx);
RRbaseSim = diff(RRtimes);

% RR interval calculation for stimulation interval
time = simOut.time(tstimIdx);
phi = simOut.Phi(tstimIdx);

% Determine indices for when phi=0 (beginning and end of heart beat)
for i = length(phi):-1:2
    if phi(i) == phi(i-1)
        phi(i-1) = [];
        time(i-1) = [];
    end
end

RRidx = find(~phi);
RRtimes = time(RRidx);
RRstimSim = diff(RRtimes);

% RR interval calculation for recovery interval
time = simOut.time(trecovIdx);
phi = simOut.Phi(trecovIdx);

% Determine indices for when phi=0 (beginning and end of heart beat)
for i = length(phi):-1:2
    if phi(i) == phi(i-1)
        phi(i-1) = [];
        time(i-1) = [];
    end
end

RRidx = find(~phi);
RRtimes = time(RRidx);
RRrecovSim = diff(RRtimes);

% Calculate percent changes in RR interval
stimDelSim = (mean(RRstimSim) - mean(RRbaseSim)) / mean(RRbaseSim)
recovDelSim = (mean(RRrecovSim) - mean(RRbaseSim)) / mean(RRbaseSim)


% Percent changes in HR
HR_vals = simOut.HR_TS;
baseHR = mean(getsampleusingtime(HR_vals,tbase(1),tbase(2)),'Weighting','time');
stimHR = mean(getsampleusingtime(HR_vals,tstim(1),tstim(2)),'Weighting','time');
recovHR = mean(getsampleusingtime(HR_vals,trecov(1),trecov(2)),'Weighting','time');
HRstimSim = (stimHR - baseHR) / baseHR
HRrecovSim = (recovHR - baseHR) / baseHR

%% bar chart plot compare percent changes in RR interval
% these values were determined by changing the stimulation inputs in the
% simulink model

% afferent only stimulation
stimAff = 0.0179;
recovAff = 0.0071;

% 1 Hz efferent and 1 Hz afferent stimulation
stimBoth = 0.1351;
recovBoth = 0.0086;

% plot formatting
fs = 18;
cyan = [0 1 1];
rose = [1 0 0.8];
red = [1 0 0];
violet = [0.5 0 1];

data = 100*[stimAff stimBoth stimDelSim stimDel; recovAff recovBoth recovDelSim  recovDel];
b = bar(data);
b(1).FaceColor = violet; %violet
b(2).FaceColor = red; %red
b(3).FaceColor = rose; %rose
b(4).FaceColor = cyan; %cyan

ylabel('RR interval % change from baseline')
legend('1 hz afferent only','1 Hz afferent, 1 Hz efferent','1 Hz afferent, 0.5 Hz efferent', 'Rajendran et al. 2019','Location','northoutside')
ylim([-5 15])
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 600, 700])
saveas(gcf,'Fig7B_RRstim_percentchange.png')


% plot RR interval traces
% tmax = trecov(2);
tplot = [0 360];
simTime         = simOut.time;
tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));

% RR interval calculation
time = simOut.time(tplotIdx);
phi = simOut.Phi(tplotIdx);

% Determine indices for when phi=0 (beginning and end of heart beat)
for i = length(phi):-1:2
    if phi(i) == phi(i-1)
        phi(i-1) = [];
        time(i-1) = [];
    end
end

RRidx = find(~phi);
RRtimes = time(RRidx);
RRinterval = diff(RRtimes);


% Plot as staircase
timePlot = repelem(RRtimes,2);
timePlot(1) = []; % get rid of first entry and last time entry to align
timePlot(end) = [];
timePlot = timePlot - tbase(1); %shift so stimulus starts at same time as data
RRintervalPlot = repelem(RRinterval,2);

%% RR interval plot
fs = 18;
lw = 2;

RRintPer = (RRint - mean(baseRR)) ./ mean(baseRR);
RRintervalPlotPer = (RRintervalPlot - mean(RRbaseSim)) ./ mean(RRbaseSim);
RRtimes = locs - 3052.2;
RRtimes(end) = [];

%% Save output for different conditions
% filename = 'AfferentOnly.mat';%'Both1Hz.mat'; %
% save(filename,'RRintervalPlotPer','timePlot')

%% plotting
% 1 hz afferent stimulus, 0.6 hz efferent
timePlot05 = timePlot;
filter_level = 0.15;
RRint05 = lowpass(RRintervalPlotPer,filter_level);

% load afferent only and both 1 Hz stimulation data and pass RR interval
% through low pass filter
load('AfferentOnly.mat')
timePlotAff = timePlot;
RRintAff = lowpass(RRintervalPlotPer,filter_level);

load('Both1Hz.mat')
timePlotBoth = timePlot;
RRintBoth = lowpass(RRintervalPlotPer,filter_level);

RRintPerFilter = lowpass(RRintPer,filter_level);

% plot
figure(3)
stairs(timePlotAff,RRintAff,'LineWidth',2,'Color',violet)
hold on
stairs(timePlotBoth,RRintBoth,'LineWidth',2,'Color',red)
stairs(timePlot05,RRint05,'LineWidth',2,'Color',rose)
stairs(RRtimes,RRintPerFilter,'LineWidth',2,'Color',cyan)
xlabel('Time (s)')
ylabel('RR interval % change from baseline')

set(gca,'FontSize',fs)
xlim([0 179])
hold off
set(gcf, 'Position',  [10, 10, 600, 500])
saveas(gcf,'Fig7C_compare_HRdataStep.png')

% save processed data in .mat file
filename = 'Aff1_Eff0_5.mat';
save(filename,'timePlot05','RRint05') % filtered RR interval data