% Script to plot figure 5 which compares ICN firing frequencies for models
% with and without a local reflex
clear; close all;

% Input parameters
Params = [2.798268 49.889025 20.086178 6.873361 0.993574 2.431832 284.973370 3.763416 1.835671 7.793137 14.090266 3.085521 14.572550 10.458438 0.183496]; %916030


%% Base model
% create simulation input object
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13';
simIn = Simulink.SimulationInput(mdlName);
simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(Params(2)), ... %output is the fmax value 
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
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(15)));

% Run simulation
simOut = sim(simIn);
%% No Local reflex
% create simulation input object
ParamsNLR = Params;
ParamsNLR(14) = 0;
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13_NLR';
simInNLR = Simulink.SimulationInput(mdlName);
simInNLR = simInNLR.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(ParamsNLR(1)), ...        
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(ParamsNLR(2)), ... 
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(ParamsNLR(3)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(ParamsNLR(4)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(ParamsNLR(5)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(ParamsNLR(6)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(ParamsNLR(7)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(ParamsNLR(8)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(ParamsNLR(9)), ... 
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(ParamsNLR(10)), ... 
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(ParamsNLR(11)), ...
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(ParamsNLR(12)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(ParamsNLR(13)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(ParamsNLR(14)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(ParamsNLR(15)));
    

% Run simulation
simOutNLR = sim(simInNLR);
%% determine indices
tlong  = [0 300];    
tmeasure = [194.8 200];
t15 = [185 200];
tplot = [150 180];

    simTime         = simOut.time;
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));
    tmeasureIdx     = find(simTime>= tmeasure(1)  & simTime <= tmeasure(2));
    t15Idx     = find(simTime>= t15(1)  & simTime <= t15(2));
    
    simTimeNLR         = simOutNLR.time;
    tplotIdxNLR        = find(simTimeNLR >= tplot(1)  & simTimeNLR <= tplot(2));
    t15IdxNLR     = find(simTimeNLR>= t15(1)  & simTimeNLR <= t15(2));

%% RR interval calculation
time = simOut.time;
phi = simOut.Phi;
Tperiod = simOut.Tperiod;


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

%% NLR RR interval
timeNLR = simOutNLR.time(tplotIdxNLR);
phiNLR = simOutNLR.Phi(tplotIdxNLR);
TperiodNLR = simOutNLR.Tperiod(tplotIdxNLR);

% Determine indices for when phi=0 (beginning and end of heart beat)
for i = length(phiNLR):-1:2
    if phiNLR(i) == phiNLR(i-1)
        phiNLR(i-1) = [];
        timeNLR(i-1) = [];
    end
end

RRidxNLR = find(~phiNLR);
RRtimesNLR = timeNLR(RRidxNLR);
RRintervalNLR = diff(RRtimesNLR);

% Plot as staircase
timePlotNLR = repelem(RRtimesNLR,2);
timePlotNLR(1) = []; % get rid of first entry and last time entry to align 
timePlotNLR(end) = [];
RRintervalPlotNLR = repelem(RRintervalNLR,2);

 
%% subpanel 1: ICN input-output firing frequencies (Healthy and MI)
fs = 18; % font size
fig = figure(2);
plotlim = [165 180]; % x-axis limits
yax = [1 12.5]; % y-axis limits


% Output firing frequencies
% PN_NA
subplot(3,2,2)
plot(simOut.time(tplotIdx),simOut.PN_NA_ff(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.PN_NA_ff(tplotIdxNLR),'r-','LineWidth',2)
title({'Output Firing Frequencies'},'FontSize',18)
xlim(plotlim)
ylim(yax)
set(gca,'FontSize',fs)

% DMV-PN firing frequency
subplot(3,2,4)
plot(simOut.time(tplotIdx),simOut.PN_DMV_ff(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.PN_DMV_ff(tplotIdxNLR),'r-','LineWidth',2)
ylim(yax)
xlim(plotlim)
set(gca,'FontSize',fs)

%LCN firing frequency
subplot(3,2,6)
plot(simOut.time(tplotIdx),simOut.LCN_ff(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.LCN_ff(tplotIdxNLR),'r-','LineWidth',2)
xlabel('Time (s)','FontSize',18)
set(gca,'FontSize',fs)
xlim(plotlim)

% Input firing frequencies
% PN_NA
subplot(3,2,1)
plot(simOut.time(tplotIdx),simOut.PN_NA_in(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.PN_NA_in(tplotIdxNLR),'r-','LineWidth',2)
ylabel({'Net'; 'activity (Hz)'},'FontSize',18)
title({'Input Firing Frequencies'},'FontSize',18)
xlim(plotlim)
ylim(yax)
set(gca,'FontSize',fs)

%PN_DMV
subplot(3,2,3)
plot(simOut.time(tplotIdx),simOut.PN_DMV_in(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.PN_DMV_in(tplotIdxNLR),'r-','LineWidth',2)
ylabel({'Net'; 'activity (Hz)'},'FontSize',18)
xlim(plotlim)
ylim(yax)
set(gca,'FontSize',fs)

% LCN
subplot(3,2,5)
plot(simOut.time(tplotIdx),simOut.LCNin(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.LCNin(tplotIdxNLR),'r-','LineWidth',2)
xlabel('Time (s)','FontSize',18)
ylabel({'Net'; 'Activity (Hz)'},'FontSize',18)
set(gca,'FontSize',fs)
xlim(plotlim)

% Format and save plot
han=axes(fig,'visible','off'); 
set(gcf, 'Position',  [10, 10, 1200, 900])
set(gca,'FontSize',fs)
saveas(gcf,'Fig5_plot_ICNff.png')

%% subpanel 2: zoomed in 
fs = 18;
plotlim = [165 180];
yax = [2.5 4.75];
figure(1)

% PN_NA output ff
subplot(2,1,1)
plot(simOut.time(tplotIdx),simOut.PN_NA_ff(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.PN_NA_ff(tplotIdxNLR),'r-','LineWidth',2)
ylabel({'Firing'; 'Frequency (Hz)'},'FontSize',18)
title({'Output Firing Frequencies'},'FontSize',18)
xlim(plotlim)
ylim(yax)
set(gca,'FontSize',fs)

% PN_DMV output ff
subplot(2,1,2)
plot(simOut.time(tplotIdx),simOut.PN_DMV_ff(tplotIdx),'b-','LineWidth',2)
hold on
plot(simOutNLR.time(tplotIdxNLR),simOutNLR.PN_DMV_ff(tplotIdxNLR),'r-','LineWidth',2)
xlabel('Time (s)','FontSize',18)
ylabel({'Firing'; 'Frequency'},'FontSize',18)
% ylim(yax)
xlim(plotlim)

% Format and save plot
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 600, 900])
saveas(gcf,'Fig5B_plot.png')
