% Script to produce figure 3B and 4B
% Plots reproduce data from Suga 1976 with input arterial pressure change
% and resulting HR and Emaxlv in open loop, and data from Greenwood 1980
% and Hainsworth 1974 (Also see Figure 2C Park 2020)


clear; close all; 

Params =  [2.798268 49.889025 20.086178 6.873361 0.993574 2.431832 284.973370 3.763416 1.835671 7.793137 14.090266 3.085521 0.3 14.572550 10.458438 0.183496];% 916030

mdlName = 'ICN_with_BR_input_model4test_clusterParams_PF2C_2lane_v2';
%% Mean arterial pressure changes
% Vary input arterial pressure in mmHg
Psa = 50:25:150;

% Preallocate
HR2lanes = zeros(length(Psa));
Emax2lanes = zeros(length(Psa));


% Simulation input with tuned model parameters and input Psa
simIn = Simulink.SimulationInput(mdlName);

% Two lanes
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
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(Params(13)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(14)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(15)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)));

% Run simulation
simOutBP = sim(simIn);

% Indices for each arterial pressure
t50 = [105 120];     % chose 15 sec range since this is what is used clinically to determine heart rate
t75 = [225 240];
t100 = [345 360];
t125 = [465 480];
t150 = [585 600];


HR2lanes(1) = mean(getsampleusingtime(simOutBP.HR,105,120),'Weighting','time');
HR2lanes(2) = mean(getsampleusingtime(simOutBP.HR,225,240),'Weighting','time');
HR2lanes(3) = mean(getsampleusingtime(simOutBP.HR,345,360),'Weighting','time');
HR2lanes(4) = mean(getsampleusingtime(simOutBP.HR,465,480),'Weighting','time');
HR2lanes(5) = mean(getsampleusingtime(simOutBP.HR,585,600),'Weighting','time');


Emax2lanes(1) = mean(getsampleusingtime(simOutBP.Emaxlv,105,120),'Weighting','time');
Emax2lanes(2) = mean(getsampleusingtime(simOutBP.Emaxlv,225,240),'Weighting','time');
Emax2lanes(3) = mean(getsampleusingtime(simOutBP.Emaxlv,345,360),'Weighting','time');
Emax2lanes(4) = mean(getsampleusingtime(simOutBP.Emaxlv,465,480),'Weighting','time');
Emax2lanes(5) = mean(getsampleusingtime(simOutBP.Emaxlv,585,600),'Weighting','time');


%% No Local reflex
Params(15) = [];
mdlName = 'ICN_with_BR_input_model4test_clusterParams_PF2C_NLR';

% Simulation input with tuned model parameters and input Psa
simIn = Simulink.SimulationInput(mdlName);

simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...        % Set ICN parameters
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
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(Params(13)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(14)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(15)));

% Run simulation
simOut = sim(simIn);


% Indices for each arterial pressure
simTime = simOut.time;

t50Idx          = find(simTime>= t50(1) & simTime <= t50(2));
t75Idx          = find(simTime>= t75(1) & simTime <= t75(2));
t100Idx         = find(simTime>= t100(1) & simTime <= t100(2));
t125Idx          = find(simTime>= t125(1) & simTime <= t125(2));
t150Idx          = find(simTime>= t150(1) & simTime <= t150(2));

HRnoLCN(1) = mean(getsampleusingtime(simOut.HR,105,120),'Weighting','time');
HRnoLCN(2) = mean(getsampleusingtime(simOut.HR,225,240),'Weighting','time');
HRnoLCN(3) = mean(getsampleusingtime(simOut.HR,345,360),'Weighting','time');
HRnoLCN(4) = mean(getsampleusingtime(simOut.HR,465,480),'Weighting','time');
HRnoLCN(5) = mean(getsampleusingtime(simOut.HR,585,600),'Weighting','time');

EmaxnoLCN(1) = mean(getsampleusingtime(simOut.Emaxlv,105,120),'Weighting','time');
EmaxnoLCN(2) = mean(getsampleusingtime(simOut.Emaxlv,225,240),'Weighting','time');
EmaxnoLCN(3) = mean(getsampleusingtime(simOut.Emaxlv,345,360),'Weighting','time');
EmaxnoLCN(4) = mean(getsampleusingtime(simOut.Emaxlv,465,480),'Weighting','time');
EmaxnoLCN(5) = mean(getsampleusingtime(simOut.Emaxlv,585,600),'Weighting','time');


%% Data from Suga 1976
HRpercent = [104.73492723492723 106.23180873180868 99.39362439362435 90.63582813582812 75.46257796257794]/ 100;
Emaxpercent = [112.18222680716232 109.95320904870681 99.39577039274926 85.84113182521554 76.9471667526343] / 100;
BaseHR = 61.2; % Used by Park 2020 based on Ursino 1998
BaseEmax = 2.695; % Used by Park 2020 based on Ursino 1998
HRSuga = BaseHR * HRpercent;
EmaxSuga = BaseEmax * Emaxpercent;



%% Plots
fig = figure(1);
ms = 6; % marker size
lw = 3; %line width
fs = 16; % font size

% Figure 3B, heart rate plot
subplot(1,2,1)
plot(Psa,HRSuga,'ko','MarkerSize',ms,'LineWidth',lw)
hold on
plot(Psa,HR2lanes(:,1),'b-','MarkerSize',ms,'LineWidth',lw)
hold on
plot(Psa,HRnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
ylabel('HR (bpm)','FontSize',fs)
ylim([25 100])
set(gca,'FontSize',fs)

%Figure 3B, elastance plot
subplot(1,2,2)
plot(Psa,EmaxSuga,'ko','MarkerSize',ms,'LineWidth',lw)
hold on
plot(Psa,Emax2lanes(:,1),'b-','MarkerSize',8,'LineWidth',lw)
hold on
plot(Psa,EmaxnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
ylabel('E_{max,lv} (mmHg/mL)','FontSize',fs)
%legend('Corresponding Experimental Data','Two Lanes','No Local Reflex','Location','northeastoutside')
ylim([1 10])
set(gca,'FontSize',fs)

% common x axis label
han = axes(fig,'visible','off');
han.XLabel.Visible = 'on';
set(gcf, 'Position',  [50, 50, 1000, 400])
saveas(gcf,'Fig3_AP.png')
hold off
%% Vary Lung Tidal Volume
clear;

mdlName = 'ICN_with_BR_input_model4test_clusterParams_ParkFig2C_Vlung_v2';
Params =  [2.798268 49.889025 20.086178 6.873361 0.993574 2.431832 284.973370 3.763416 1.835671 7.793137 14.090266 3.085521 0.3 14.572550 10.458438 0.183496];% 916030

Vlung = [0, 0.75, 1.25, 1.75]; % Percent increase from nominal

% Park model params
BRclosedLoop_model_initialize_params;

% Preallocate
HR2lanes = zeros(length(Vlung));
Emax2lanes = zeros(length(Vlung));


% Simulation input with tuned model parameters and input Psa
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
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(Params(13)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(14)), ... 
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(Params(15)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(16)));

% Run simulation
simOut = sim(simIn);

% Indices for each arterial pressure
t50 = [105 120];     % chose 15 sec range since this is what is used clinically to determine BPM
t75 = [225 240];
t100 = [345 360];
t125 = [465 480];


HR2lanes(1) = mean(getsampleusingtime(simOut.HR,105,120),'Weighting','time');
HR2lanes(2) = mean(getsampleusingtime(simOut.HR,225,240),'Weighting','time');
HR2lanes(3) = mean(getsampleusingtime(simOut.HR,345,360),'Weighting','time');
HR2lanes(4) = mean(getsampleusingtime(simOut.HR,465,480),'Weighting','time');


Emax2lanes(1) = mean(getsampleusingtime(simOut.Emaxlv,105,120),'Weighting','time');
Emax2lanes(2) = mean(getsampleusingtime(simOut.Emaxlv,225,240),'Weighting','time');
Emax2lanes(3) = mean(getsampleusingtime(simOut.Emaxlv,345,360),'Weighting','time');
Emax2lanes(4) = mean(getsampleusingtime(simOut.Emaxlv,465,480),'Weighting','time');

%% % No Local reflex
Params(15) = [];
mdlName = 'ICN_with_BR_input_model4test_clusterParams_ParkFig2C_VlungNLR';

% Simulation input with tuned model parameters 
simIn = Simulink.SimulationInput(mdlName);

simIn = simIn.setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(Params(1)), ...        % Set ICN parameters
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
    [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(Params(13)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(Params(14)), ...
    [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(15)));

% Run simulation
simOut = sim(simIn);

HRnoLCN(1) = mean(getsampleusingtime(simOut.HR,105,120),'Weighting','time');
HRnoLCN(2) = mean(getsampleusingtime(simOut.HR,225,240),'Weighting','time');
HRnoLCN(3) = mean(getsampleusingtime(simOut.HR,345,360),'Weighting','time');
HRnoLCN(4) = mean(getsampleusingtime(simOut.HR,465,480),'Weighting','time');


EmaxnoLCN(1) = mean(getsampleusingtime(simOut.Emaxlv,105,120),'Weighting','time');
EmaxnoLCN(2) = mean(getsampleusingtime(simOut.Emaxlv,225,240),'Weighting','time');
EmaxnoLCN(3) = mean(getsampleusingtime(simOut.Emaxlv,345,360),'Weighting','time');
EmaxnoLCN(4) = mean(getsampleusingtime(simOut.Emaxlv,465,480),'Weighting','time');

%% Experimental data
HRGreenwood = [68.67647058823529 77.6470588235294  82.94117647058822 90.14705882352939];
HRPark = [76.02847548690396 80.21383478844861 84.92975151108124 93.75795836131628];
xGreenwood = [0, 0.73015873015872930, 1.25396825396825480, 1.78571428571428650];
xEmax = [0, 0.75, 1.25, 1.75];
EmaxPark = [2.3947848045674136 2.2945652173913045 2.3663043478260875 2.378072024593764];
EmaxGreenwood = [ 2.2539174352217834 2.0771739130434783 2.1228238910847606 2.173722002635046];
%% Plots
ms = 6;
lw = 3;
fs = 16;

fig = figure(2);

% Figure 4B, heart rate plot
subplot(1,2,1)
plot(xGreenwood,HRGreenwood,'ko','MarkerSize',ms,'LineWidth',lw)
hold on
plot(xGreenwood,HRnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
hold on
plot(Vlung,HR2lanes(:,1),'b-','MarkerSize',ms,'LineWidth',lw)
ylabel('HR (bpm)','FontSize',fs)
ylim([30 100])
ax = gca;
set(gca,'FontSize',fs)

% Figure 4B, elastance plot
subplot(1,2,2)
plot(xEmax,EmaxGreenwood,'ko','MarkerSize',ms,'LineWidth',lw)
hold on
plot(xEmax,EmaxnoLCN,'r-','MarkerSize',ms,'LineWidth',lw)
hold on
plot(Vlung,Emax2lanes(:,1),'b-','MarkerSize',ms,'LineWidth',lw)
ylabel('E_{max,lv} (mmHg/mL)','FontSize',fs)
ylim([1 10])
set(gca,'FontSize',fs)


% common x axis label
han = axes(fig,'visible','off');
han.XLabel.Visible = 'on';
set(gcf, 'Position',  [50, 50, 1000, 400])
hold off
saveas(gcf,'Fig4_Vlung.png')