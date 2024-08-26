function Mastitskaya_plot_bar_HRMAP(preIR_MAP_models, postIR_MAP_models, preIR_HR_models, postIR_HR_models,altMdlName)
% preIR_MAP_exp: Pre-I/R mean arterial pressure data for experimental data (1x1 vector)
% postIR_MAP_exp: Post-I/R mean arterial pressure data for experimental data (1x1 vector)
% preIR_HR_exp: Pre-I/R heart rate data for experimental data (1x1 vector)
% postIR_HR_exp: Post-I/R heart rate data for experimental data (1x1 vector)
% preIR_MAP_models: Pre-I/R mean arterial pressure data for models (1xn vector)
% postIR_MAP_models: Post-I/R mean arterial pressure data for models (1xn vector)
% preIR_HR_models: Pre-I/R heart rate data for models (1xn vector)
% postIR_HR_models: Post-I/R heart rate data for models (1xn vector)

%% Experimental data
% Data from Mastitskaya 2012
GeeHealthy = [69 91]; % model baseline HR and MAP
ratControl_PI = [408 106]; %Mastitskaya control pre-ischemia HR and MAP

% Mastitskaya percent changes from control pre-ischemia to experimental
% group pre-ischemia
ratHR_PI = [408 401 420 379];
ratMAP_PI = [106 115 99 99];
HRchange_control_PI = (ratHR_PI - ratControl_PI(1)) / ratControl_PI(1);
MAPchange_control_PI = (ratMAP_PI - ratControl_PI(2)) / ratControl_PI(2);

% Calculate scale factor and scaled HR and MAP
scale = ratControl_PI ./ GeeHealthy;

scaledMastitskayaHR = ratControl_PI(1) ./ scale(1);
scaledMastitskayaMAP = ratControl_PI(2) ./ scale(2);

% Calculate experimental group pre-ischemia values based on percent change
preIR_HR_exp = scaledMastitskayaHR * (1 + HRchange_control_PI);
preIR_MAP_exp = scaledMastitskayaMAP * (1 + MAPchange_control_PI);

% Calculate end-ischemia HR and MAP based on percent change
% Data
EI_HR = [405 387 393 384]; % 120 min reperfusion % end-ischemia[411 421 415 392];
EI_MAP = [91 82 83 71]; % 120 min reperfusion % end-ischemia [101 104 79 76];
% Percent change
HRchange_PI_EI = (EI_HR - ratHR_PI) ./ ratHR_PI;
MAPchange_PI_EI = (EI_MAP - ratMAP_PI) ./ ratMAP_PI;
% Final value
postIR_HR_exp = preIR_HR_exp .* (1+HRchange_PI_EI);
postIR_MAP_exp = preIR_MAP_exp .* (1+MAPchange_PI_EI);

% Final vectors
% HRconverted = [HRconverted_PI; HRconverted_EI]; % control, DMV silencing, DMV activation, DMV activation + atropine
% MAPconverted = [MAPconverted_PI; MAPconverted_EI];

% Error bars
% raw values
HR_SEM = [8 9 6 14; 9 13 9 10]; % 120 min reperfusion [preischemia;endischemia] % end-ischemia [8 9 6 14; 9 13 8 15];
MAP_SEM = [7 9 6 10; 4 7 4 3]; % 120 min reperfusion % end-ischemia[7 9 6 10; 7 10 4 6];

% Convert standard error to 95% confidence intervals
numAnimals = 7;
HR_SD = HR_SEM .* sqrt(numAnimals); % convert to standard deviation
MAP_SD = MAP_SEM .* sqrt(numAnimals);
HR_95CI = HR_SD .* 1.96;
MAP_95CI = MAP_SD .* 1.96;

% convert to percent
HRerr_per = HR_95CI ./ [ratHR_PI; EI_HR];
MAPerr_per = MAP_95CI./ [ratMAP_PI; EI_MAP];

% final error bar values converted
% HRerr_converted = HRerr_per .* HRconverted;
% MAPerr_converted = MAPerr_per .* MAPconverted;

%preIR MAP 69-130
preIR_MAP_exp = 99.5;
std_preIR_MAP_exp = 30.5; % 95% CI used %MAPerr_per(1,1) * preIR_MAP_exp; % 25.5;%
std_postIR_MAP_exp = MAPerr_per(2,1) * postIR_MAP_exp;
std_preIR_HR_exp = HRerr_per(1,1) * preIR_HR_exp;
std_postIR_HR_exp = HRerr_per(2,1) * postIR_HR_exp;

% convert to RR interval in ms
preIR_RR_exp = 940;% move to mideel fo error bars since HR conversion makes it look wonky 1000*60/65.5;%.*60./preIR_HR_exp; %
std_preIR_RR_exp = 1000*60/10.5;%HRerr_per(1,1) .* preIR_RR_exp; % 55-76 bpm -->
% postIR range 58-116 BPM
postIR_RR_exp = 780;% move to middle of error bars since HR conversion makes it look wonky 1000.*60./postIR_HR_exp;%1000*60/(58+19);
std_postIR_RR_exp = HRerr_per(2,1) .* postIR_RR_exp;%1000*60/19;%

% Bounds based on Osculati for postIR
postIR_RR_exp = 780;% move to middle of error bars since HR conversion makes it look wonky %1000*60/87; %690 ms
std_postIR_RR_exp = 250;%1000*60/29; %58-116 bpm --> 0.53-1.03 s/beat
postIR_MAP_exp = 86.9;
std_postIR_MAP_exp = 17.6; % 95% CI not stddev

%% Create figure with 1x2 subplot grid
figure;
fs = 16;
rng default
max_jitter = 0.1;
min_jitter = 0;
x_jitter = (max_jitter - min_jitter).*randn(length(preIR_MAP_models),1) + min_jitter;
disp('xjitter')
size(x_jitter)
color = [0 0 1 0.5];
lineColor = [0    0.4471    0.7412];
barColor = [0.7608    0.3961    0.3961];

% Plot experimental data for MAP
subplot(1, 2, 1);
hold on;
bar(1, preIR_MAP_exp(1),'FaceColor',barColor)
bar(2, postIR_MAP_exp(1), 'FaceColor',barColor,'HandleVisibility','off');
% errorbar(x,y,neg,pos)
errorbar([1, 2], [preIR_MAP_exp(1), postIR_MAP_exp(1)], [std_preIR_MAP_exp(1), std_postIR_MAP_exp(1)], 'k', 'LineStyle', 'none','HandleVisibility','off','LineWidth',6);
ylabel('Mean arterial pressure (mm Hg)');
xtick_label = {'Pre-I/R'; 'Post-I/R'};
set(gca,'FontSize',fs,'XTick',[1:2],'xticklabel',xtick_label)
xlim([0.5, 2.5]);
ylim([40 145])

% Plot experimental data for heart rate
subplot(1, 2, 2);
hold on;
bar(1, preIR_RR_exp(1), 'FaceColor',barColor,'HandleVisibility','off');
bar(2, postIR_RR_exp(1), 'FaceColor',barColor,'HandleVisibility','off');
errorbar([1, 2], [preIR_RR_exp(1), postIR_RR_exp(1)], [150, 250], [150, 250],'k', 'LineStyle', 'none','HandleVisibility','off','LineWidth',6);
ylabel('RR interval (ms)');
set(gca,'FontSize',fs,'XTick',[1:2],'xticklabel',xtick_label)
xlim([0.5, 2.5]);
ylim([520 1200])

ms = 4;
% Plot model data for MAP
subplot(1, 2, 1);
hold on;
for i = 1:length(preIR_MAP_models)
    plot([1+x_jitter(i), 2+x_jitter(i)], [preIR_MAP_models(i), postIR_MAP_models(i)], 'o-', 'Color', lineColor,'MarkerSize',ms,'MarkerFaceColor','b','MarkerEdgeColor','k');
end
ttl = title('C');
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.27; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  


% Plot model data for heart rate
preIR_RR_models = 1000.*60./preIR_HR_models;
postIR_RR_models = 1000.*60./postIR_HR_models;
subplot(1, 2, 2);
hold on;
for i = 1:length(preIR_HR_models)
    plot([1+x_jitter(i), 2+x_jitter(i)], [preIR_RR_models(i), postIR_RR_models(i)], 'o-', 'Color', lineColor,'MarkerSize',ms,'MarkerFaceColor','b','MarkerEdgeColor','k');
end
ttl = title('D');
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.27; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

% Legend
subplot(1, 2, 1);
legend('Experimental', 'Model');


set(gcf, 'Position',  [10, 10, 800, 400])
filename = ['./plot_HRMAP_062924' altMdlName '.png']; % /plots
saveas(gcf,filename)


end




