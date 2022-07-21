% Script to plot sensitivity analysis results in the order fmin, fmax,
% fmid, k for PN_NA, PN_DMV, and LCN. The input gains are the final three
% parameters and in the order k_fev, k_fesh, k_BR
clear; close all;
Xlab = {'1','2','3'};


% colors 
orange = [1 0.5 0]; % for map
purple = [0.4940 0.1840 0.5560];
lightBlue = [0.3010 0.7450 0.9330];

%Load HR data
load('PAWN_HR_916030.mat')
M = 15;
MM = 1:1:M;


% sort by into user-specified order
[tmp,ranking] = sort(-PAWN_median(end,:));

% ranking = 1:1:M; % add/remove this line to put in consistent order
% order NA, DMV, LCN, LCN input gains
% min, max, mid, k
% fevEmax, fesh, BR gain
ranking = [ 1 2 3 4 9 10 11 12 5 6 7 8 13 15 14];

hfig= figure; 
fs = 24 ;
b = bar(PAWN_median(end,ranking),'FaceColor',[0.3010 0.7450 0.9330]) ;


axis([0,M+1,0,0.6])
hold on
plot([1:M],ones(1,M)*KS_dummy_mean(end),'r','LineWidth',3) % plots red line
%text(3,0.4,['N = ' int2str(N_red(end)) ', n = ' num2str(n) ],'FontSize',fs)
set(gca,'FontSize',fs)
for i=1:M; Xlab{i}=num2str(ranking(i)); end
set(gca,'XTick',[1:M],'XTickLabel',Xlab,'YTick',[0,0.2,0.4],'FontSize',16)
%text(1.45,0.5,['Mean Arterial Pressure'],'FontSize',fs) % Title
%for i=1:M; text(i-0.5,PAWN_median_ub(end,ranking(i))+0.02,num2str(ranking(i)),'FontSize',20); end
ylabel('Sensitivity Index','FontSize',fs)
xlabel('Parameters','FontSize',fs)
for i = 1:M % connect upper and lower bound with a line
    plot([i i],[PAWN_median_lb(end,ranking(i)) PAWN_median_ub(end,ranking(i))],'k','LineWidth',1.5)
end
set(hfig, 'Position', [0 0 1600 400])
saveas(gcf,'PAWNbar_HR_916030.png')
