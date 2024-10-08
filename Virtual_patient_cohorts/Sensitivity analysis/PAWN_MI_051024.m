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

% Purpose: identify parameters that significantly affect HR of those used
% to model myocardial ischemia for Mastitskaya 2012 simulations
% Parameter values of interest: 
%  {'PN_{NA} f_{min}','PN_{NA} f_{max}','PN_{NA} f_{mid}','PN_{NA} k', ...
%     'LCN f_{min}','LCN f_{max}','LCN f_{mid}','LCN k', ...
%     'PN_{DMV} f_{min}','PN_{DMV} f_{max}','PN_{DMV} f_{mid}','PN_{DMV} k', '\tau_{PNDMV}', ...
%     'k_{fevEmax}','k_{BR}','k_{CP}','k_{fesh}', ...
%     'BR f_{min}','BR f_{max}','BR f_{mid}','BR k', ...
%     'CP f_{min}','CP f_{max}','CP f_{mid}','CP k', ...
%     'LSR f_{min}','LSR f_{max}','LSR f_{mid}','LSR k', ...
%     'NA f_{min}','NA f_{max}','NA f_{mid}','NA k', ...
%     'NActr f_{min}','NActr f_{max}','NActr f_{mid}','NActr k', ...
%     'DMV f_{min}','DMV f_{max}','DMV f_{mid}','DMV k', ...
%     'C_{la}',V_{ula}','R_{la}','P_{0lv}','k_{Elv}','V_{ulv}','k_{Rlv}','E_{maxlv0}','fes_{inf}','k_a',
%     'T_{sys,0}', 'k_{sys}', 'G_{T,v}', 'G_{T,s}'};
% These values were chosen based on
% the ones that Park et al. 2020 retuned to model heart failure

% Results: for MAP and HR parameters 5, 9, 4, 1 are significant. THese are
% Cla, P0lv, kElv, fes_inf
clear; close all;
%% Step 1: Add paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/PAWN/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/')
addpath(genpath('C:\Users\mmgee\MATLAB\safe_R1.1'))

% date for filenames
date = '20240510';
%% Step 2: setup the cardiovascular control model

% Define uncertain inputs (parameters):
% ICN parameters
PNDMVdelay = 0.3; 
ICNparams = [1.6906 7.3295 13.811104 3.230162 1.988567 18.375719 521.217197 3.003231 2.415242 17.718990 14.183806 13.356069 PNDMVdelay 3.329861 2.661685 5.642977 0.066794];

% NA/DMV
NAparams = [4.88, 15.78, 59.83, 23, 0.61, 11, 12.81, 7, 2.5901, 6.66, 42.91, 33.5]; %  NA, NActr, DMV

% NTS
NTSparams = [0.30, 21.50, 37.07, 21, 0.45, 28.33, 10.2, 7, 2.75, 31.57, 11.13, 2]; %  BR, CPR, LSR; fmin, fmax, fmid, k

% Cardiac parameters
Cla = 19.23;
Vula = 25;
Rla = 2.5e-3;
kRlv = 3.75e-4;
P0lv = 1.5;
kElv = 0.014;
Vulv = 16.77;
Emaxlv0 = 1.283;
fes_inf = 2.1;

% Baroreceptor
ka = 11.758;

% Systolic duration
Tsys0 = 0.4; 
ksys = 0.075; 

% Gain of sympathetic and parasympathetic inputs to heart rate
Gts = 0.13;
Gtv = 0.09;


MIparams = [ICNparams(1), ICNparams(2), ICNparams(3), ICNparams(4), ... %fmin_NA_PN, fmax_NA_PN, midptNA_PN, kNA,...
    ICNparams(5), ICNparams(6), ICNparams(7), ICNparams(8), ... %fmin_LCN, fmax_LCN, midptLCN,kLCN, ...
    ICNparams(9), ICNparams(10), ICNparams(11), ICNparams(12), ICNparams(13) ... %fmin_DMV_PN, fmax_DMV_PN, midpt_DMV_PN, kDMV_PN, tau_DMV_PN, ...
    ICNparams(14), ICNparams(15), ICNparams(16), ICNparams(17), ... %LCN_fevEmaxgain, LCN_BRgain, LCN_CPgain, LCN_feshgain, ...
    NTSparams(1), NTSparams(2), NTSparams(3), NTSparams(4),... %fmin_NTS_BR, fmax_NTS_BR, midpt_NTS_BR, k_NTS_BR, ...
    NTSparams(5), NTSparams(6), NTSparams(7), NTSparams(8), ... %fmin_NTS_CP, fmax_NTS_CP, midpt_NTS_CP, k_NTS_CP, ...
    NTSparams(9), NTSparams(10), NTSparams(11), NTSparams(12), ... %fmin_NTS_LS, fmax_NTS_LS, midpt_NTS_LS, k_NTS_LS, ...
    NAparams(1), NAparams(2), NAparams(3), NAparams(4), ... %fmin_NA, fmax_NA, midpt_NA, k_NA, ...
    NAparams(5), NAparams(6), NAparams(7), NAparams(8), ... %fmin_NActr, fmax_NActr, midpt_NActr, k_NActr, ...
    NAparams(9), NAparams(10), NAparams(11), NAparams(12), ... %fmin_DMV, fmax_DMV, midpt_DMV, k_DMV, ...
    Cla, Vula, Rla, P0lv, kElv, Vulv, kRlv, Emaxlv0, fes_inf, ...
    ka, ...
    Tsys0, ksys, Gts, Gtv];

M = length(MIparams); % number of inputs
labelparams= {'PN_{NA} f_{min}','PN_{NA} f_{max}','PN_{NA} f_{mid}','PN_{NA} k', ...
    'LCN f_{min}','LCN f_{max}','LCN f_{mid}','LCN k', ...
    'PN_{DMV} f_{min}','PN_{DMV} f_{max}','PN_{DMV} f_{mid}','PN_{DMV} k', '\tau_{PNDMV}', ...
    'k_{fevEmax}','k_{BR}','k_{CP}','k_{fesh}', ...
    'BR f_{min}','BR f_{max}','BR f_{mid}','BR k', ...
    'CP f_{min}','CP f_{max}','CP f_{mid}','CP k', ...
    'LSR f_{min}','LSR f_{max}','LSR f_{mid}','LSR k', ...
    'NA f_{min}','NA f_{max}','NA f_{mid}','NA k', ...
    'NActr f_{min}','NActr f_{max}','NActr f_{mid}','NActr k', ...
    'DMV f_{min}','DMV f_{max}','DMV f_{mid}','DMV k', ...
    'C_{la}','V_{ula}','R_{la}','P_{0lv}','k_{Elv}','V_{ulv}','k_{Rlv}','E_{maxlv0}','fes_{inf}','k_a', ...
    'T_{sys,0}', 'k_{sys}', 'G_{T,v}', 'G_{T,s}'};

% Baseline model parameters
kRSA = 0.5;
BPswitch = -1; % closed loop


% Vary parameter values over a 2 fold range from nominal. 
factor = 1.5;
xmin = MIparams./factor;
xmax = MIparams.*factor;


distrpar=cell(M,1); for i=1:M; distrpar{i}=[xmin(i) xmax(i)]; end

%% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);


%% Define name of simulink model
mdlName     = 'ICN_model_v15_Mastitskaya2012_control_r2020b';

%% Step 3: Apply PAWN
% Number of samples and distribution function for sampling
N = 45000; n = 6; %2100
N_MAP = N;
DistrFun = 'unif' ;


X = AAT_sampling('lhs',M,DistrFun,distrpar,N);
X_MAP = X;



% Feed X into model to get model output of interest Yu (CO, HR, MAP, etc)
% Create an array of simulation input objects and specify the sweep value for each simulation
simIn(1:N) = Simulink.SimulationInput(mdlName);
for idx = 1:N
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on');%, ...
%         'TimeOut', 240);

simIn(1,idx) = simIn(1,idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(X(idx,1)), ...
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
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/kRSA'], 'Value', num2str(kRSA), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmin'], 'Value', num2str(X(idx,18)), ...         % NTSparams
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmax'], 'Value', num2str(X(idx,19)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/fmid'], 'Value', num2str(X(idx,20)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/BR cell group/k, gain'], 'Value', num2str(X(idx,21)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmin'], 'Value', num2str(X(idx,22)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmax'], 'Value', num2str(X(idx,23)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/fmid'], 'Value', num2str(X(idx,24)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/CP cell group/k, gain'], 'Value', num2str(X(idx,25)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmin'], 'Value', num2str(X(idx,26)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmax'], 'Value', num2str(X(idx,27)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/fmid'], 'Value', num2str(X(idx,28)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NTS/LS cell group/k, gain'], 'Value', num2str(X(idx,29)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmin'], 'Value', num2str(X(idx,30)), ...     % NAparams
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmax'], 'Value', num2str(X(idx,31)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/fmid'], 'Value', num2str(X(idx,32)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA/k, gain'], 'Value', num2str(X(idx,33)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmin'], 'Value', num2str(X(idx,34)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmax'], 'Value', num2str(X(idx,35)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/fmid'], 'Value', num2str(X(idx,36)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/NA-ctr/k, gain'], 'Value', num2str(X(idx,37)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmin'], 'Value', num2str(X(idx,38)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmax'], 'Value', num2str(X(idx,39)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/fmid'], 'Value', num2str(X(idx,40)), ...
        [mdlName '/Autonomic Nervous System/Brainstem Vagal efferents/DMV/k, gain'], 'Value', num2str(X(idx,41)), ...
        [mdlName '/Left Heart/Cla'], 'Value', num2str(X(idx,42)), ...   
        [mdlName '/Autonomic Nervous System/Vula'], 'Value', num2str(X(idx,43)), ...
        [mdlName '/Left Heart/Rla'], 'Value', num2str(X(idx,44)), ...% heart_params2
        [mdlName '/Autonomic Nervous System/Left ventricle/P0lv'], 'Value', num2str(X(idx,45)), ...
        [mdlName '/Autonomic Nervous System/Left ventricle/kElv'], 'Value', num2str(X(idx,46)), ...
        [mdlName '/Autonomic Nervous System/Left ventricle/Vulv'], 'Value', num2str(X(idx,47)), ...
        [mdlName '/Left Heart/Subsystem2/kRlv'], 'Value', num2str(X(idx,48)), ...
        [mdlName '/Autonomic Nervous System/Emaxlv_inv/Emaxlv0_denervated (Emaxlv0)'], 'Value', num2str(X(idx,49)), ...
        [mdlName '/Autonomic Nervous System/Symp. Efferent Pathways/fes_inf'], 'Value', num2str(X(idx,50)), ...
        [mdlName '/Autonomic Nervous System/Carotid Sinus (baroreceptors)/ka'], 'Value', num2str(X(idx,51)), ...
        [mdlName '/Autonomic Nervous System/Heart_period/Phi_activity2/Tsys0'], 'Value', num2str(X(idx,52)), ...
        [mdlName '/Autonomic Nervous System/Heart_period/Phi_activity2/ksys'], 'Value', num2str(X(idx,53)), ...
        [mdlName '/Autonomic Nervous System/Heart_period/Gv'], 'Value', num2str(X(idx,54)), ...
        [mdlName '/Autonomic Nervous System/Heart_period/Gts'], 'Value', num2str(X(idx,55)), ...
        [mdlName '/Autonomic Nervous System/BPswitch'], 'Value', num2str(BPswitch));
end

%% Run simulations
simout = parsim(simIn);
% simout = sim(simIn(1,1));

%% Analysis for HR
meanHR = zeros(N,1); % Preallocate HR output vector
MAP = zeros(N,1); % Preallocate HR output vector
for i = 1:N
    % Determine indices to use for 15 seconds at steady state
    simTime         = simout(1,i).time;
    tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
    physOutputs = PhysOutputs_Gen_TS(simout(1,i),tshort,tshort);
    SP = physOutputs(1);
    DP = physOutputs(2);
    % Calculate avg HR over 15 seconds
    HR = mean(getsampleusingtime(simout(1,i).HR_TS,tshort(1),tshort(2)));
    if length(HR) == 0 % remove outputs for simulations that failed
        meanHR(i,1) = NaN;
    else
        meanHR(i,1) = HR;
        MAP(i,1) = (SP + 2*DP)/3;
    end
end

%% HR
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
% set(gca,'FontSize',fs,'XTick',round(xc2{i},1),'YLim',[0,1],'XGrid','on')
xlabel('x_1'); ylabel('KS'); 

filename = ['cdfandKS_HR_' date '.png'];
saveas(gcf,filename)

%% Barchart of KS stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PAWN sensitivity indices for different choices of N
% and different aggregation statistics
% (this will take some computing time!!!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 10 ;
N = N - numRemovedVals;
N_red = [  2000 3000 4000 N ] ; % size of reduced sample for approximation of PAWN indices

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
set(hfig, 'Position', [0 0 3600 400])
filename = ['newPAWNbar_HR_' date '.png'];
saveas(gcf,filename)


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
set(hfig, 'Position', [0 0 3600 400])
filename = ['newPAWNbarCompare_HR_' date '.png'];
saveas(gcf,filename)

% Save output in .mat file
filename = ['PAWN_HR_' date '.mat'];
save(filename,'PAWN_median','KS_dummy_mean','PAWN_median_lb','PAWN_median_ub')


%% MAP
% reset variables used in HR PAWN analysis
N = N_MAP;
X = X_MAP;
Y = MAP; % vector (N,1)

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
% set(gca,'FontSize',fs,'XTick',round(xc2{i},1),'YLim',[0,1],'XGrid','on')
xlabel('x_1'); ylabel('KS'); 

filename = ['cdfandKS_MAP_' date '.png'];
saveas(gcf,filename)

%% Barchart of KS stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PAWN sensitivity indices for different choices of N
% and different aggregation statistics
% (this will take some computing time!!!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 10 ;
N = N - numRemovedVals;
N_red = [ 2000 3000 4000 N  ] ;  % size of reduced sample for approximation of PAWN indices

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
filename = ['newPAWNbar_MAP_' date '.png'];
saveas(gcf,filename)


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
set(hfig, 'Position', [0 0 3600 400])
filename = ['newPAWNbarCompare_MAP_' date '.png'];
saveas(gcf,filename)

% Save output in .mat file
filename = ['PAWN_MAP_' date '.mat'];
save(filename,'PAWN_median','KS_dummy_mean','PAWN_median_lb','PAWN_median_ub')

%% Exit code
close_system
delete(myPool);
exit