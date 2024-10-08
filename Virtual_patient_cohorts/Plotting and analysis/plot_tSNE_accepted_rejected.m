% tSNE of parameter values for post-IR models colored by accepted versus
% rejected parameter sets
% Just try it with postIR baroreceptor model for now
% Inputs: 
%   paramMatrix: matrix of 28 model parameters for accepted and rejected
%   postIR models
%   labels: labels are either accepted or rejected
clear; close all;
my_dir = pwd;
addpath(genpath(my_dir))
rng default
% Load postIR baroreceptor data
altMdlName = 'all';
load(['postIR_' altMdlName '_06232024.mat'])

% Load postIR baroreceptors accepted/rejected indices
load(['postIRacceptedIdx_' altMdlName '.mat'],'indices','data')

% Format paramMatrix from postIRstore
numPatients = 59;
paramMatrix = zeros(numPatients*100,51);
labels = repmat({'Rejected'},numPatients*100,1);

for i = 1%1:numPatients
    % Get the parameter matrix for the current element
    idx = (i-1)*100 + 1;
    i = 9;
    paramMatrix(idx:idx+99,:) = postIRstore(i).params;

    % Label accepted versus rejected
    % Accepted: rows where columns 1, 2, and 6 are 1
    % These criteria are SP, HR, and Osculati baroreflex curve
    crit = postIRstore(i).crit; % 100x8 array for patient i
    acceptedRows = find(crit(:, 1) == 1 & crit(:, 2) == 1 & crit(:, 6) == 1);

    % Update the labels for accepted rows
    if ~isempty(acceptedRows)
        labels(idx + acceptedRows - 1) = {'Accepted'}; % Adjust index based on accepted rows
    end

end
%% Plot tsne
perplexity = 5;
Y = tsne(paramMatrix(1:100,:),'Algorithm','barneshut','Distance','cityblock','Perplexity',perplexity,'Standardize',1); 
figure;
fs = 20;
clr = {'b' 'r'};
sym = 'oo';
siz = 6;
% Separate indices for Accepted and Rejected labels
acceptedIdx = strcmp(labels(1:100), 'Accepted');
rejectedIdx = strcmp(labels(1:100), 'Rejected');

% Plot the points
hold on;
scatter(Y(rejectedIdx, 1), Y(rejectedIdx, 2), siz, clr{2}, 'filled');  % Rejected (red, filled circles)
scatter(Y(acceptedIdx, 1), Y(acceptedIdx, 2), siz, clr{1}, 'filled');  % Accepted (blue, filled circles)
hold off;

xlabel('t-SNE 1','FontSize',fs)
ylabel('t-SNE 2','FontSize',fs)
legend('Rejected','Accepted','Location', 'best'); 
set(gca,'FontSize',fs)
filename = ['tSNE_06232024_' altMdlName '_acceptedRejected.png'];
saveas(gcf,filename)

%% Plot PCA
xNames = {'PN_NA_fmin','PN_NA_fmax','PN_NA_fmid','PN_NA_k','LCN_fmin','LCN_fmid','LCN_k',...
    'PN_DMV_fmin','PN_DMV_fmid','LCN_CPgain','NTS_BR_fmid','NTS_CP_fmax',...
    'NTS_LS_fmax','NTS_LS_fmid','NA_fmin','NA_fmax','NA_fmid','NA_ctr_fmid','DMV_fmax','DMV_fmid',...
    'C_{la}','P_{0lv}','k_{Elv}','V_{ulv}','k_{Rlv}','E_{maxlv0}','fes_{inf}','k_a'};
XSNames = labels;
paramMatrix_norm = zscore(paramMatrix(1:100,[1 2 3 4 5 7 8 9 11 16 20 23 27 28 30 31 32 36 39 40 42 45 46 47 48 49 50 51]),0,1);
[coeff,score,latent,tsquared,explained]  = pca(paramMatrix_norm);
figure;
b = biplot(coeff(:,1:2),'Scores',score(:,1:2),...
    'varlabels',xNames,...
    'ObsLabels',XSNames(901:1000),'MarkerSize',16);

% Change obs color
bID = get(b,'tag');
bPt = b(strcmp(bID,'obsmarker')); 
for i = 1:numel(bPt)
    switch XSNames{i}
        case {'Rejected'}
            set(bPt(i),'Color',[1 0 0 0.5]); % red
        case {'Accepted'}
            set(bPt(i),'Color',[0 1 0]); % blue
    end
end

% legend('Rejected','Accepted','Location', 'best'); 
set(gca,'FontSize',fs)
set(gcf,'Position',[5 5 700 600])
filename = ['pca_06232024_' altMdlName '_acceptedRejected.png'];
saveas(gcf,filename)
