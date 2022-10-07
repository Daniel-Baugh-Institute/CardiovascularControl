% Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
% Michelle Gee
% October 3, 2022

% Script to produce Figure 3

clear; close all; 
% Input parameters
Params = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 3.329861 2.661685 5.642977 0.066794];
kRSA = 0.5;

% Define simulink input
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

% Run simulation
simOut = sim(simIn);

% Calculate cardiovascular metrics
disp('Systolic BP, diastolic BP, HR, CO, EF')
PhysOutputs_Gen_TS(simOut,[165 180], [150 180])
%% Plot
plotlim = [165 180];
simTime         = simOut.time; % time output vector from simulation
tplotIdx       = find(simTime>= plotlim(1) & simTime <= plotlim(2));

% Plot formatting
fs = 16; % font size
lw = 2; % line width
rows = 2; % subplot dimensions
cols = 1;

% PV loop
fig = figure(1);
plot(simOut.Vlv(tplotIdx),simOut.Plv(tplotIdx),'b','LineWidth',lw)
ylabel({'Left ventricular'; 'pressure (mm Hg)'},'FontSize',fs)
xlabel('Left ventricular volume (mL)','FontSize',fs)
ax = gca;
set(gca,'FontSize',fs)
xlim([0 150])
ylim([0 150])
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 400, 400])
saveas(gcf,'Fig3_PV.png')

% Mean arterial pressure and lung tidal volume
figure(2)
subplot(rows,cols,1)
plot(getsampleusingtime(simOut.Psa_TS,plotlim(1),plotlim(2)),'b','LineWidth',lw)
ylabel({'Mean arterial'; 'pressure (mm Hg)'},'FontSize',fs)
title('')
xlabel('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)
ylim([70 130])

subplot(rows,cols,2)
plot(getsampleusingtime(simOut.Vlung_TS,plotlim(1),plotlim(2)),'b','LineWidth',lw)
ylabel({'Lung tidal'; 'volume (L)'},'FontSize',fs)
xlabel('Time (s)','FontSize',fs)
title('')
ax = gca;
set(gca,'FontSize',fs)
xlim(plotlim)


% save figure
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 400, 600])
saveas(gcf,'Fig3_psa_Vlung.png')