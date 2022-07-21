clear; close all;
%% load data and format into array
load('RSA_GPR_NLR_100_HR10916030_v13.mat')
var_NLR = var_base;
load('RSA_GPR_noRSA_100_HR10916030_v13.mat')
var_NRSA = var_base;
load('RSA_GPR_NLR_NRSA_100_HR10916030_v13.mat')
var_NLR_NRSA = var_base;
load('RSA_GPR_base_100_HR10916030_v13.mat')

var_all = [var_base' var_NLR' var_NRSA' var_NLR_NRSA'];
%% Analysis
% remove NAN values if necesary
var_all_compare = var_all;
var_all(any(isnan(var_all), 2), :) = [];

% calculate information entropy
var_all = log(var_all);
var_all = -var_all;
mdl_mean = mean(var_all,1);

% T-test to determine if mean likelihood is different across model
% structures. Zero means not significantly different, one means
% significantly different
h12 = ttest(var_all(:,1),var_all(:,2))
h13 = ttest(var_all(:,1),var_all(:,3))
h14 = ttest(var_all(:,1),var_all(:,4))
h23 = ttest(var_all(:,2),var_all(:,3))
h24 = ttest(var_all(:,2),var_all(:,4))
h34 = ttest(var_all(:,3),var_all(:,4))

% select 5 examples to plot
var_ex = [var_all(12:13,:); var_all(43:45,:)];

% plot
fs = 14;
xvalues = {'Local reflex and RSA gate','No local reflex','No RSA gate','No local reflex, no RSA gate'};
yvalues = 1:1:length(var_ex(:,1))';
h = heatmap(xvalues,yvalues,var_ex,'Ylabel','ICN parameter set','Colormap',parula,'CellLabelColor','none');
set(gca,'FontSize',fs)
set(gcf, 'Position',  [10, 10, 800, 500])
saveas(gcf,'plot_RSA_heatmap.png')