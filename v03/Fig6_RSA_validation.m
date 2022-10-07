% Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
% Michelle Gee
% October 3, 2022

% Script to produce Figure 6

% Script to determine RSA amplitude for inputs of lung tidal volumes of
% 0.5-1.6 L in increments of 0.1 L. Runs files
% ICN_with_BR_input_model4test_clusterParams_ICNtune_v15.slx
% sequentially to calculate RSA amplitude

% calculates HR from RR interval as done by Kobayashi rather than directly
% calculating RSA amplitude from instantaneous HR

% Compares model results to Kobayashi 1998 

%% Run model and calculate RSA amplitude
clear; close all;

%vector of lung tidal volumes
Vlung = 0.4:0.1:1.6;
alpha = 0.08:0.02:0.4;
Vln = 2.22:0.08:3.5;

% Preallocate vectors
RSAamplitude = zeros(1,length(Vlung));
stdRSA = zeros(1,length(Vlung));
minRSA = zeros(1,length(Vlung));
maxRSA = zeros(1,length(Vlung));

% simulation input parameters
mdlName = 'ICN_model_v15';

for i = 1:length(Vlung)
    Params = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 3.329861 2.661685 5.642977 0.066794];
   kRSA = 0.5;
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
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(kRSA), ...
        [mdlName '/Subsystem/alpha'], 'Value', num2str(alpha(i)), ...
        [mdlName '/Subsystem/Vln'], 'Value', num2str(Vln(i)));
    
     disp('Simulating lung tidal volume of:')
    Vlung(i)
    
    % Run simulation
   simOut(1,i) = sim(simIn);
    
    % determine indices for 5.2 second measurement windows based on
    % Kobayashi 1998 measurement windows
    tlong  = [0 300];     % plotting timerange
    tmeasure = [194.8 200];
    window = 5.2;
    tstart = 185:window:181+window*length(Vlung);
    tend = tstart + 5.2;
    
    % RR interval calculation
    simTime         = simOut(1,i).time;
    tplotIdx        = find(simTime>= tlong(1)  & simTime <= tlong(2));
    
    time = simTime(tplotIdx);
        phi = simOut(1,i).Phi(tplotIdx);
        
        % Determine indices for when phi=0 (beginning and end of heart beat)
        for k = length(phi):-1:2
            if phi(k) == phi(k-1)
                phi(k-1) = [];
                time(k-1) = [];
            end
        end
        
        RRidx = find(~phi);
        RRtimes = time(RRidx);
        RRint = diff(RRtimes);
        RRint(1) = [];
        RRint(end) = []; % remove last RR interval value because it is usually an incomplete measurement
        HR = 60./RRint;
        
    % calculate RR interval for each measurement window
    for j = 1:7
        tmeasureIdx     = find(simTime>= tstart(j)  & simTime <= tend(j));
        HRIdx = find(RRtimes >= tstart(j)  & RRtimes <= tend(j));
        HRwindow = HR(HRIdx);
        RSA(1,j) = max(HRwindow) - min(HRwindow);
      
    end

    stdRSA(i) = std(RSA);
    minRSA(i) = min(RSA);
    maxRSA(i) = max(RSA);
    RSAamplitude(i) = mean(RSA);
end

RSAamplitude
stdRSA

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


%% Plot
% plot formatting
fs = 18; % font size
lw = 2; % line width
ms = 9; % marker size

% Error bars
err = stdRSA;

figure(22)
h1 = errorbar(Vlung,RSAamplitude,err,'bo','MarkerFaceColor','b','MarkerSize',ms); % simulation
hold on
plot(Vlung,RSAmdl,'b--','LineWidth',lw) % simulation line of best fit
h2 = errorbar(RSAtidalK,RSAampK,exp_sd,'co','MarkerFaceColor','c','MarkerSize',ms); %kobayashi 1998
plot(RSAtidalK,RSAk,'c--','LineWidth',lw) % kobayashi line of best fit
xlabel('Lung tidal volume (L)','FontSize',18)
ylabel('RSA amplitude (bpm)','FontSize',18)
xlim([0.38 1.65])

legend([h1,h2],'Model prediction, R^{2}=0.93','Kobayashi et al. 1998, R^{2}=0.98','Location','northoutside')%,'Freyschuss experimental data (0.1 Hz), R^{2}=0.99','Ben-tal model (0.20 Hz), R^{2}=0.99') % ,'Freyschuss experimental data (0.1 Hz)','Freyschuss experimental fit'
hold off

set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 600, 600])
saveas(gcf,'Fig6_RSA_validation.png')
