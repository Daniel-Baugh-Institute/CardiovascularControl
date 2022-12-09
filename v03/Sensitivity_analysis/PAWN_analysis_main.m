% PAWN global sensitivity analysis for cardiovascular control model

% This script was adapted from the PAWN sensitivity analysis approach

% Original code can be found: https://www.safetoolbox.info/pawn-method/
% The PAWN toolbox must be downloaded to run this file.

% (Pianosi and Wagener, 2018) Pianosi F, Wagener T. Distribution-based
% sensitivity analysis from a generic input-output sample. Environmental
% Modelling & Software. 2018;108:197-207. 

% Pianosi, F. and Wagener, T. (2015), A simple and efficient method for
% global sensitivity analysis based on cumulative distribution functions,
% Env. Mod. & Soft., 67, 1-11.


clear; close all;
%% Step 1: Add paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2: setup the cardiovascular control model

% Define uncertain inputs (parameters):
M = 18 ; % number of inputs
labelparams= { '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15', '16', '17', '18'  } ; % input names
%{ 'minval_NA_PN','fmax_NA_PN','midptNA_PN','kNA','minval_LCN','fmax_LCN','midptLCN','kLCN','minval_DMV_PN','fmax_DMV_PN','midptDMV_PN','kDMV_PN','tau_DMV_PN','LCN_fevEmaxgain','LCN_BRgain', 'LCN_CPgain','LCN_feshgain','RSA_gain'  } ; % input names

% Baseline model parameters
PNDMVdelay = 0.3;
Params = [0.857242 12.624676 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 PNDMVdelay 3.329861 2.661685 5.642977 0.066794];
kRSA = 0.5;

% Vary parameter values over a 3 fold range from nominal. kRSA (Params(18))
% was varied from 0 to 1 because it represents the relative decrease in
% input to the NA during inhalation which represents gating of the inputs
% to the NA observed in respiratory sinus arrhythmia.
xmin = Params./3;
xmin(18) = 0;
xmax = Params.*3;
xmax(18) = 1;

distrpar=cell(M,1); for i=1:M; distrpar{i}=[xmin(i) xmax(i)]; end

%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);


%% Define name of simulink model
mdlName     = 'ICN_model_v15';

%% Step 3: Apply PAWN
% Number of samples and distribution function for sampling
N = 3300; n = 6;
DistrFun = 'unif' ;


X = AAT_sampling('lhs',M,DistrFun,distrpar,N);


% Feed X into model to get model output of interest Yu (CO, HR, MAP, etc)
% Create an array of simulation input objects and specify the sweep value for each simulation
simIn(1:length(X)) = Simulink.SimulationInput(mdlName);
for idx = 1:length(X)
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);

simIn(1,idx) = simIn(1,idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(X(idx,1)), ...        % Set ICN parameters
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(X(idx,2)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(X(idx,3)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(X(idx,4)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(X(idx,5)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(X(idx,6)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(X(idx,7)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(X(idx,8)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(X(idx,9)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(X(idx,10)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(X(idx,11)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(X(idx,12)), ...
        [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(X(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(X(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(X(idx,15)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_CPgain'], 'Value', num2str(X(idx,16)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(X(idx,17)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(X(idx,18)));
end

% Run simulations
simout = parsim(simIn);

% Analysis for HR
meanHR = zeros(length(X),1); % Preallocate HR output vector
for i = 1:length(X)
    % Determine indices to use for 15 seconds at steady state
    simTime         = simout(1,i).time;
    tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
    
    % Calculate avg HR over 15 seconds
    HR = mean(getsampleusingtime(simout(1,i).HR_TS,tshort(1),tshort(2)));
    if length(HR) == 0 % remove outputs for simulations that failed
        meanHR(i,1) = NaN;
    else
        meanHR(i,1) = HR;
    end
end


Y = meanHR; % vector (N,1)

% Check Y for NaN and count number of NaN simulations
z = isnan(Y);
k = find(z); % indices of NaN values of Y

% Number of NaN values reomoved
numRemovedVals = length(k);
sprintf('Number of NaN simulation is %1.0f',numRemovedVals)

% Remove NaN values from X and Y
for i = length(k):-1:1
    Y(k(i)) = [];
    X(k(i),:) = [];
end


[KS,KS_dummy,YF,Fu,Fc,YY2,xc2,NC,XX2,YU,idx_bootstrap] = pawn_ks_givendata(X,Y,n);

%% Plotting cdf and KS
col=gray(n+1); % Color for colormap, (n+1) so that the last line is not white
map = colormap(col(1:n,:));% Remove white color
i =1; % input factor under investigation
fs = 22; ms = 10 ; % font size and marker size

% (g): plot CDFs
subplot(1,2,1)
for k=1:n
    plot(YF,Fc{i,k},'Color',map(k,:),'LineWidth',3.5); hold on
end
plot(YF,Fu,'r','LineWidth',3.5)
set(gcf, 'Position',  [100, 100, 1200, 500])
set(gca,'FontSize',fs); xlabel('y'); ylabel('cdf'); 
% (h): plot KS
subplot(1,2,2)
for k=1:n
    plot(xc2{i}(k),KS(k,i),'ok','MarkerFaceColor',map(k,:),'MarkerSize',ms); hold on
end
set(gca,'FontSize',fs,'XTick',round(xc2{i},1),'YLim',[0,1],'XGrid','on')
xlabel('x_1'); ylabel('KS'); 

saveas(gcf,'cdfandKS.png')

%% Barchart of KS stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PAWN sensitivity indices for different choices of N
% and different aggregation statistics
% (this will take some computing time!!!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 10 ;
N = N - numRemovedVals;
N_red = [ 2200 2700 2800 3000 ] ; % size of reduced sample for approximation of PAWN indices

nboot = 50  ;
alfa  = 0.05 ;

PAWN_median   = nan(length(N_red),M); PAWN_median_lb = PAWN_median; PAWN_median_ub = PAWN_median;
PAWN_mean     = nan(length(N_red),M); PAWN_mean_lb = PAWN_mean; PAWN_mean_ub = PAWN_mean;
PAWN_max      = nan(length(N_red),M); PAWN_max_lb = PAWN_max; PAWN_max_ub = PAWN_max;
KS_dummy_mean = nan(1,length(N_red));

for r=1:length(N_red) % WARNING: running this for loop may take some time!
    
    idx = randperm(N,N_red(r)) ;
    [KS_median,KS_mean,KS_max,KS_dummy] = pawn_indices_givendata(X(idx,:),Y(idx),n,nboot) ;
    
    % Take statistics across bootstrap resamples:
    PAWN_median(r,:) = mean(KS_median) ;
    median_lb = sort(KS_median) ; median_lb = median_lb(max(1,round(nboot*alfa/2)),:);
    median_ub = sort(KS_median) ; median_ub = median_ub(round(nboot*(1-alfa/2)),:);
    PAWN_median_lb(r,:) = median_lb ;
    PAWN_median_ub(r,:) = median_ub ;
    %
    PAWN_mean(r,:)   = mean(KS_mean)  ;
    mean_lb = sort(KS_mean) ; mean_lb = mean_lb(max(1,round(nboot*alfa/2)),:);
    mean_ub = sort(KS_mean) ; mean_ub = mean_ub(round(nboot*(1-alfa/2)),:);
    PAWN_mean_lb(r,:) = mean_lb ;
    PAWN_mean_ub(r,:) = mean_ub ;
    %
    PAWN_max(r,:)   = mean(KS_max)  ;
    max_lb = sort(KS_max) ; max_lb = max_lb(max(1,round(nboot*alfa/2)),:);
    max_ub = sort(KS_max) ; max_ub = max_ub(round(nboot*(1-alfa/2)),:);
    PAWN_max_lb(r,:) = max_lb ;
    PAWN_max_ub(r,:) = max_ub ;
    %
    KS_dummy_mean(r) = mean(KS_dummy)  ;
    %
    fprintf('iteration %d (N_red=%d) of %d completed \n',r,N_red(r),length(N_red))

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 4: Bar plot of PAWN indices for given (N,n) [using stat=median] 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
[tmp,ranking] = sort(-PAWN_median(end,:));
%X_Labels_swat_ranked=X_Labels_swat(ranking);


hfig= figure; 
fs = 24 ;
b = bar(PAWN_median(end,ranking)) ;
set(b,'FaceColor',[126 126 126]/256)
axis([0,M+1,0,0.6])
hold on
plot([1:M],ones(1,M)*KS_dummy_mean(end),'r','LineWidth',2)
text(3,0.4,['N = ' int2str(N_red(end)) ', n = ' num2str(n) ],'FontSize',fs)
set(gca,'FontSize',fs)
for i=1:M; Xlab{i}=num2str(ranking(i)); end
set(gca,'XTick',[1:M],'XTickLabel',Xlab,'YTick',[0,0.2,0.4],'FontSize',16)
%for i=1:M; text(i-0.5,PAWN_median_ub(end,ranking(i))+0.02,num2str(ranking(i)),'FontSize',20); end
ylabel('PAWN index','FontSize',fs)
xlabel('input factors','FontSize',fs)
for i = 1:M % connect upper and lower bound with a line
    plot([i i],[PAWN_median_lb(end,ranking(i)) PAWN_median_ub(end,ranking(i))],'k','LineWidth',1.5)
end
set(hfig, 'Position', [0 0 1600 400])
saveas(gcf,'newPAWNbar_HR_supp.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: Bar plots of PAWN indices for given n and varying N
% [using stat=median or mean or max ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hfig= figure; fs = 24 ;
stat = 1 ; % median (creates Fig. 6)
%stat = 2 ; % mean (creates Fig. 8)
%stat = 3 ; % max (creates Fig. 9)

for r=1:length(N_red)
    
    %subaxis(length(N_red),1,r, 'SpacingVert', 0.01,'SpacingHoriz',0.01,'MarginRight',0.01,'MarginTop',0.01,'MarginBottom',0.05);
    subplot(length(N_red),1,r)
    
    if stat==1; b = bar(PAWN_median(r,ranking)) ; ylabel('PAWN index','FontSize',fs) ; end%ylabel('median(KS)','FontSize',fs); end
    if stat==2; b = bar(PAWN_mean(r,ranking))   ; ylabel('mean(KS)','FontSize',fs)  ; end
    if stat==3; b = bar(PAWN_max(r,ranking))    ; ylabel('max(KS)','FontSize',fs)   ; end
    set(b,'FaceColor',[126 126 126]/256)
    axis([0,M+1,0,0.6])
    hold on
    plot([1:M],ones(1,M)*KS_dummy_mean(r),'r','LineWidth',2)
    text(3,0.4,['N = ' int2str(N_red(r)) ', n = ' num2str(n) ],'FontSize',fs)
    
    % add axis ticks and label
    for i=1:M; Xlab{i}=num2str(ranking(i)); end
    set(gca,'FontSize',fs,'YTick',[0,0.2,0.4])
    if r<length(N_red) 
        set(gca,'XTick',[1:M],'XTickLabel',{}); 
    else
        set(gca,'XTick',[1:M],'XTickLabel',Xlab,'FontSize',16)
        xlabel('input factors','FontSize',fs);
    end
    
    for i = 1:M % connect upper and lower bound with a line
        if stat==1; plot([i i],[PAWN_median_lb(r,ranking(i)) PAWN_median_ub(r,ranking(i))],'k','LineWidth',1.5); end
        if stat==2; plot([i i],[PAWN_mean_lb(r,ranking(i)) PAWN_mean_ub(r,ranking(i))],'k','LineWidth',1.5);     end
        if stat==3; plot([i i],[PAWN_max_lb(r,ranking(i)) PAWN_max_ub(r,ranking(i))],'k','LineWidth',1.5);       end        
    end
    
end
set(hfig, 'Position', [0 0 1500 1600])
saveas(gcf,'newPAWNbarCompare_HR_supp.png')

% Save output in .mat file
save('PAWN_HR_v15.mat','PAWN_median','KS_dummy_mean','PAWN_median_lb','PAWN_median_ub')





%% Exit code
close_system
delete(myPool);
exit