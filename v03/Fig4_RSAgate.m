% Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
% Michelle Gee
% October 3, 2022

% Script to produce Figure 4
clear; close all; 

% Input parameters
Params = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 3.329861 2.661685 5.642977 0.066794];
kRSA = 0.5;
%% Define simulink input
mdlName = 'ICN_model_v15';
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


% analysis
tmax = 197;
tplot = [118 tmax];
simTime         = simOut.time; % time output vector from simulation
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
RRintervalPlot = repelem(RRinterval,2);

%% No RSA model
% define simulation input
mdlName = 'ICN_model_v15';
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
    [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(1)); 
    % this sets kRSA = 1, effectively removing the RSA gate
    
% run simulation
simOutNRSA = sim(simIn);


% analysis
tmax = 197;
tplot = [118 tmax];
simTime         = simOutNRSA.time;
tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));

% RR interval calculation
time = simOutNRSA.time(tplotIdx);
phi = simOutNRSA.Phi(tplotIdx);

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
timePlotNRSA = repelem(RRtimes,2);
timePlotNRSA(1) = []; % get rid of first entry and last time entry to align 
timePlotNRSA(end) = [];
RRintervalPlotNRSA = repelem(RRinterval,2);


%% Plot
% Plot formatting
plotlim = [164 179];
fs = 12; % font size
lw = 2; % line width
rows = 4; % subplot dimensions
cols = 2;

fig = figure(1);
% RSA gate, NA activity
subplot(rows,cols,1)
plot(getsampleusingtime(simOut.fevHR,plotlim(1),plotlim(2)),'LineWidth',lw)
ylabel({'NA net'; 'activity (Hz)'},'FontSize',fs)
title('')
xlabel('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

% No RSA gate, NA activity
subplot(rows,cols,2)
plot(getsampleusingtime(simOutNRSA.fevHR,plotlim(1),plotlim(2)),'LineWidth',lw)
ylabel({''},'FontSize',fs)
title('')
xlabel('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

% RSA gate, PN_NA activity
subplot(rows,cols,3)
plot(getsampleusingtime(simOut.PN_NA_ff,plotlim(1),plotlim(2)),'LineWidth',lw)
ylabel({'PN_{NA} net'; 'activity (Hz)'},'FontSize',fs)
xlabel('')
title('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

% no RSA gate, PN_NA activity
subplot(rows,cols,4)
plot(getsampleusingtime(simOutNRSA.PN_NA_ff,plotlim(1),plotlim(2)),'LineWidth',lw)
ylabel('','FontSize',fs)
xlabel('')
title('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

% RSA gate, RR interval
subplot(rows,cols,5)
plot(timePlot,RRintervalPlot,'LineWidth',2)
ylabel('RR interval (s)','FontSize',fs)
xlabel('Time (s)','FontSize',fs)
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)
ylim([0.8 0.95])

% no RSA gate, RR interval
subplot(rows,cols,6)
plot(timePlotNRSA,RRintervalPlotNRSA,'LineWidth',2)
ylabel('','FontSize',fs)
xlabel('Time (s)','FontSize',fs)
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)
ylim([0.8 0.95])

% RSA gate, lung tidal volume (not shown in figure but used to align
% inspiration and exhalation intervals)
subplot(rows,cols,7)
plot(getsampleusingtime(simOut.Vlung_TS,plotlim(1),plotlim(2)),'LineWidth',lw)
ylabel({'Lung tidal'; 'volume (L)'},'FontSize',fs)
xlabel('Time (s)','FontSize',fs)
title('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

% no RSA gate, lung tidal volume (not shown in figure but used to align
% inspiration and exhalation intervals)
subplot(rows,cols,8)
plot(getsampleusingtime(simOutNRSA.Vlung_TS,plotlim(1),plotlim(2)),'LineWidth',lw)
ylabel({''},'FontSize',fs)
xlabel('Time (s)','FontSize',fs)
title('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)

% save figure
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 800, 800])
saveas(gcf,'Fig4_RSAgate.png')

