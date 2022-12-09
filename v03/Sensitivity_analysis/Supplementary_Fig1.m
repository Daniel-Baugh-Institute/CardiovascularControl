%% Code to plot Supplementary Figure 1, the results of PAWN sensitivity analysis
% Michelle Gee
% December 9, 2022

clear; close all;
Xlab = {'1','2','3'};
% Load HR data
load('PAWN_HR_v15.mat')

% Number of parameter values
M = 18;
MM = 1:1:M;


% Reorder parameters for preferred order in the plot
[tmp,ranking] = sort(-PAWN_median(end,:));

% ranking = 1:1:M; % add/remove this line to put in consistent order
% order NA, DMV, LCN, LCN input gains, kRSA
% min, max, mid, k
% fBR, fCP, fevEmax, fesh, RSA gain
ranking = [ 1 2 3 4 9 10 11 12 13 5 6 7 8 15 16 14 17 18];

hfig= figure; 
fs = 24 ; % font size

% Plot
b = bar(PAWN_median(end,ranking),'FaceColor',[0.3010 0.7450 0.9330]) ;

% Format location on axis
axis([0,M+1,0,0.6])
hold on
plot([1:M],ones(1,M)*KS_dummy_mean(end),'r','LineWidth',3) % plots red line

% Plot labeling and formatting 
set(gca,'FontSize',fs)
for i=1:M; Xlab{i}=num2str(ranking(i)); end
set(gca,'XTick',[1:M],'XTickLabel',Xlab,'YTick',[0,0.2,0.4],'FontSize',16)
ylabel('Sensitivity Index','FontSize',fs)
xlabel('Parameters','FontSize',fs)

for i = 1:M % connect upper and lower bound with a line
    plot([i i],[PAWN_median_lb(end,ranking(i)) PAWN_median_ub(end,ranking(i))],'k','LineWidth',1.5)
end

% Set figure size and position
set(hfig, 'Position', [0 0 1600 400])

% Save figure
saveas(gcf,'Supplementary_Fig1.png')


