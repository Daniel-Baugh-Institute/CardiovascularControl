% Script to determine RSA amplitude for inputs of lung tidal volumes of
% 0.5-2 L in increments of 0.1 L. Runs files
% ICN_with_BR_input_model4test_clusterParams_ICNtune_v13.slx
% sequentially to calculate RSA amplitude

% Compares model results to Kobayashi 1998 and Ben-Tal 2012 and
% Freyschuss 1976 data 

clear; close all;

%% Run model and calculate RSA amplitude
%vector of lung tidal volumes
Vlung = 0.5:0.1:2;
alpha = 0.1:0.02:0.4;
Vln = 2.3:0.08:3.5;



mdlName = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13_RSA_25Hz';

% Preallocate vectors
RSAamplitude = zeros(1,length(Vlung));

% Input parameters
for i = 1:length(Vlung)
    Params = [2.798268 49.889025 20.086178 6.873361 0.993574 2.431832 284.973370 3.763416 1.835671 7.793137 14.090266 3.085521 14.572550 10.458438 0.183496];
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
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(Params(15)), ...
        [mdlName '/Subsystem/alpha'], 'Value', num2str(alpha(i)), ...
        [mdlName '/Subsystem/Vln'], 'Value', num2str(Vln(i)));
    
 
    % Run simulation
    simOut(1,i) = sim(simIn);
    
    % determine indices
    tlong  = [0 300];     % plotting timerange
    tmeasure = [194.8 200];
    tstart = 152:6:194;
    tend = tstart + 5.2;
    
    % RSA Amplitude calculation
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tlong(1)  & simTime <= tlong(2));
    
    % measure RSA amplitude at different time points
    for j = 1:length(tstart)
        tmeasureIdx     = find(simTime>= tstart(j)  & simTime <= tend(j));
        RSA(j) = max(simOut(1,i).HR(tmeasureIdx)) - min(simOut(1,i).HR(tmeasureIdx));
        
    end
    stdRSA(i) = std(RSA);
    minRSA(i) = min(RSA);
    maxRSA(i) = max(RSA);
    RSAamplitude(i) = mean(RSA);
end

RSAAmplitudenoLR = RSAamplitude 

% Linear regression
mdl = fitlm(Vlung,RSAamplitude)
T = mdl.Coefficients;
coeffs = T{1:2,{'Estimate'}};
RSAmdl = coeffs(2)*Vlung + coeffs(1);


%% Plotting Comparison to other models and experimental data
% Kobayashi data
% Lung tidal volume values
RSAtidalK = [0.391837
0.497856
0.600187
0.69521
0.797595
0.894399
0.991316
1.097305
1.194261
1.298348
1.393366
1.49228
1.587215];

% RSA amplitude
RSAampK = [1.572581
2.110215
2.849462
3.521505
3.891129
4.899194
5.134409
5.873656
5.840054
7.083333
7.788978
6.88172
8.158602];

% standard deviation Kobayashi
exp_sd = [0.840054
1.612903
2.049731
2.520161
2.419355
2.520162
2.788979
2.688172
2.822581
2.688172
3.59543
2.184139
3.091398];

% Comparison to Kobayashi fit
RSAk = 5.45*RSAtidalK + 0.31;


% Freyschuss
VlungF = [1 1.5 2]; % Lung volumes
RSAampF = [11.9 16.9 22.6]; % RSA amplitude
sd_F = [1.9 1.7 2]; % standard error of the mean
% line of best fit from excel: y = 10.7x + 1.0833
RSAfitF = 10.7*VlungF + 1.0833;

% Ben-Tal
VlungBT = [0.33657917
0.713378493
1.113463167
1.515664691
1.919983065
2.324301439]; % Lung tidal volumes

RSAampBT = [11.70460026
13.3137547
15.30236123
17.51245994
20.00734199
22.72374302]; %RSA amplitude

% Comparison to Ben-tal fit
RSAbt = 5.55*VlungBT + 9.44; % Fit from Ben-Tal


% Plot
% Error bars
err = stdRSA;
fs = 18; % font size
lw = 1; % line width
medgreen = [0 0.5 0]; % color for plotting

figure(22)
h1 = errorbar(Vlung,RSAamplitude,err,'bo','MarkerFaceColor','b');
hold on
plot(Vlung,RSAmdl,'b--','LineWidth',lw)
h2 = errorbar(RSAtidalK,RSAampK,exp_sd,'co','MarkerFaceColor','c');
plot(RSAtidalK,RSAk,'c--','LineWidth',lw)
h3 = errorbar(VlungF,RSAampF,sd_F,'o','Color',medgreen,'MarkerFaceColor',medgreen);
plot(VlungF,RSAfitF,'--','Color',medgreen,'LineWidth',lw)
h4 = plot(VlungBT,RSAampBT,'mo','MarkerFaceColor','m');
plot(VlungBT,RSAbt,'m--')

xlabel('Lung Tidal Volume (L)','FontSize',18)
ylabel('RSA Amplitude (bpm)','FontSize',18)
legend([h1,h2,h3,h4],'Model prediction (0.25 Hz), R^{2}=0.91','Kobayashi experimental data (0.25 Hz), R^{2}=0.98','Freyschuss experimental data (0.1 Hz), R^{2}=0.99','Ben-tal model (0.20 Hz), R^{2}=0.99','Location','northeastoutside') % ,'Freyschuss experimental data (0.1 Hz)','Freyschuss experimental fit'
hold off
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 1200, 600])
saveas(gcf,'Fig2_plot_RSA.png')
