% PAWN global sensitivity analysis for cardiovascular control model

% This script provides an application example
% of the PAWN sensitivity analysis approach (Pianosi and Wagener, 2015)
%
% MODEL AND STUDY AREA
%
%
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method
% for global sensitivity analysis based on cumulative distribution
% functions, Env. Mod. & Soft., 67, 1-11.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk
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
M = 15 ; % number of inputs
labelparams= { '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'};%,'16'  } ; % input names
%{ 'minval_NA_PN','fmax_NA_PN','midptNA_PN','kNA','minval_LCN','fmax_LCN','midptLCN','kLCN','minval_DMV_PN','fmax_DMV_PN','midptDMV_PN','kDMV_PN','tau_DMV_PN','LCN_fevEmaxgain','LCN_BRgain','LCN_feshgain'  } ; % input names
Params = [2.798268 49.889025 20.086178 6.873361 0.993574 2.431832 284.973370 3.763416 1.835671 7.793137 14.090266 3.085521 14.572550 10.458438 0.183496];
xmin = Params./2;
xmax = Params.*2;

distrpar=cell(M,1); for i=1:M; distrpar{i}=[xmin(i) xmax(i)]; end

% Load model output:
% Set up parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
myCluster.JobStorageLocation = getenv('TMPDIR');
myPool = parpool(myCluster, myCluster.NumWorkers);


% Define name of simulink model
mdlName     = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v13';
%% Step 3: Apply PAWN
% (edit) these values
N = 5000; n = 6;
DistrFun = 'unif' ;


X = AAT_sampling('lhs',M,DistrFun,distrpar,N);


% Feed X into model to get model output of interest Yu (CO, HR, MAP, etc)
% Create an array of simulation input objects and specify the sweep value for each simulation
simIn(1:length(X)) = Simulink.SimulationInput(mdlName);
for idx = 1:length(X)
    simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
        'SaveOutput', 'on', ...
        'TimeOut', 240);

    % Two lanes
%     simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(X(idx,1)), ...        % Set ICN parameters
%         [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(X(idx,2)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(X(idx,3)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(X(idx,4)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(X(idx,5)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(X(idx,6)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(X(idx,7)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(X(idx,8)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(X(idx,9)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(X(idx,10)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(X(idx,11)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(X(idx,12)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(X(idx,13)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(X(idx,14)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(X(idx,15)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(X(idx,16)));
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
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(X(idx,13)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(X(idx,14)), ...
        [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(X(idx,15)));
end

% Run simulations
simout = parsim(simIn);
%simout = sim(simIn(1,1));

% Analysis for HR
meanHR = zeros(length(X),1);
for i = 1:length(X)
    % Determine indices to use for 15 seconds at steady state
    simTime         = simout(1,i).time;
    tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
    tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
    
    % Calculate avg HR over 15 seconds
    meanHR(i,1)         = mean(simout(1,i).HR(tshortIdx));
end


Y = meanHR; % vector (2000,1)

% Analysis for MAP
% meanPsa = zeros(length(X),1);
% for i = 1:length(X)
%     % Determine indices to use for 15 seconds at steady state
%     simTime         = simout(1,i).time;
%     tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
%     tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
%     
%     % Calculate avg HR over 15 seconds
%     meanPsa(i,1)         = mean(simout(1,i).Psa(tshortIdx));
% end
% 
% 
% Y = meanPsa; % vector (2000,1)

% Check Y for NaN and count number of NaN simulations
z = isnan(Y);
k = find(z); % indices of NaN values of Y

% Number of Nan values reomoved
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
N_red = [ 2000 3000 4000 4500 ] ; % size of reduced sample for approximation of PAWN indices

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
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 4: Bar plot of PAWN indices for given (N,n) [using stat=median] 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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
saveas(gcf,'newPAWNbar_HR_916030.png')

%%

% %% Step 2: setup the cardiovascular control model
% 
% % Define uncertain inputs (parameters):
% M = 16 ; % number of inputs
% labelparams= { '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15', '16'  } ; % input names
% %{ 'minval_NA_PN','fmax_NA_PN','midptNA_PN','kNA','minval_LCN','fmax_LCN','midptLCN','kLCN','minval_DMV_PN','fmax_DMV_PN','midptDMV_PN','kDMV_PN','tau_DMV_PN','LCN_fevEmaxgain','LCN_BRgain','LCN_feshgain'  } ; % input names
% 
% % Base
% % parameter ranges: from 2x and 1/2x the values for the actual parameters: 2.351891 17.267573 22.389733 9.615043 1.433934 30.442230 142.647356 9.451185 2.030215 29.134814 22.745358 3.289855 3860.427321 4.143256 1.729508 0.09510
% %2x xmin = [1.17594550000000,8.63378650000000,11.1948665000000,4.80752150000000,0.716967000000000,15.2211150000000,71.3236780000000,4.72559250000000,1.01510750000000,14.5674070000000,11.3726790000000,1.64492750000000,1930.21366050000,2.07162800000000,0.864754000000000,0.0475500000000000];
% % 3x min[7.05567300000000,51.8027190000000,67.1691990000000,28.8451290000000,4.30180200000000,91.3266900000000,427.942068000000,28.3535550000000,6.09064500000000,87.4044420000000,68.2360740000000,9.86956500000000,11581.2819630000,12.4297680000000,5.18852400000000,0.285309000000000];
% % 2x xmax = [4.70378200000000,34.5351460000000,44.7794660000000,19.2300860000000,2.86786800000000,60.8844600000000,285.294712000000,18.9023700000000,4.06043000000000,58.2696280000000,45.4907160000000,6.57971000000000,7720.85464200000,8.28651200000000,3.45901600000000,0.190200000000000];
% % 3x max[0.783963666666667,5.75585766666667,7.46324433333333,3.20501433333333,0.477978000000000,10.1474100000000,47.5491186666667,3.15039500000000,0.676738333333333,9.71160466666667,7.58178600000000,1.09661833333333,1286.80910700000,1.38108533333333,0.576502666666667,0.0317010000000000];
% xmax = [3.52783650000000,25.9013595000000,33.5845995000000,14.4225645000000,2.15090100000000,45.6633450000000,213.971034000000,14.1767775000000,3.04532250000000,43.7022210000000,34.1180370000000,4.93478250000000,5790.64098150000,6.21488400000000,2.59426200000000,0.142650000000000]; %1.5
% xmin = [1.56792733333333,11.5117153333333,14.9264886666667,6.41002866666667,0.955956000000000,20.2948200000000,95.0982373333333,6.30079000000000,1.35347666666667,19.4232093333333,15.1635720000000,2.19323666666667,2573.61821400000,2.76217066666667,1.15300533333333,0.0634000000000000]; %1.5
% 
% % Homog ICN
% % Homogeneous ICN parameters [0.595743 16.493073 41.171402 1.625914 5.611725 17.130081 146.692651 6.551255 4.671538 26.576990 25.318147 7.094506 6505.901731 6.577870 18.739766 0.394420];
% %xmin = [0.397162000000000,10.9953820000000,27.4476013333333,1.08394266666667,3.74115000000000,11.4200540000000,97.7951006666667,4.36750333333333]; % 1.5x
% %xmax = [0.893614500000000,24.7396095000000,61.7571030000000,2.43887100000000,8.41758750000000,25.6951215000000,220.038976500000,9.82688250000000]; %1.5x
% % xmax = [0.744678750000000,20.6163412500000,51.4642525000000,2.03239250000000,7.01465625000000,21.4126012500000,183.365813750000,8.18906875000000,5.83942250000000,33.2212375000000,31.6476837500000,8.86813250000000,8132.37716375000,8.22233750000000,23.4247075000000,0.493025000000000];
% % xmin = [0.476594400000000,13.1944584000000,32.9371216000000,1.30073120000000,4.48938000000000,13.7040648000000,117.354120800000,5.24100400000000,3.73723040000000,21.2615920000000,20.2545176000000,5.67560480000000,5204.72138480000,5.26229600000000,14.9918128000000,0.315536000000000];
% 
% % No local Circuit parameters [2.37467000000000,22.0421310000000,23.7918360000000,6.04086500000000,0.841829000000000,44.5346680000000,121.110560000000,5.44541500000000,2.23984600000000,25.8539820000000,58.0509740000000,7.53819700000000,6041.70122200000,9.92038100000000,0.495271000000000]; % noLCN;
% %xmin = [1.18733500000000,11.0210655000000,11.8959180000000,3.02043250000000,0.420914500000000,22.2673340000000,60.5552800000000,2.72270750000000,1.11992300000000,12.9269910000000,29.0254870000000,3.76909850000000,3020.85061100000,4.96019050000000,0.247635500000000];%2x
% %xmax = [4.74934000000000,44.0842620000000,47.5836720000000,12.0817300000000,1.68365800000000,89.0693360000000,242.221120000000,10.8908300000000,4.47969200000000,51.7079640000000,116.101948000000,15.0763940000000,12083.4024440000,19.8407620000000,0.990542000000000]; %2x
% % xmin = [1.58311333333333,14.6947540000000,15.8612240000000,4.02724333333333,0.561219333333333,29.6897786666667,80.7403733333333,3.63027666666667,1.49323066666667,17.2359880000000,38.7006493333333,5.02546466666667,4027.80081466667,6.61358733333333,0.330180666666667];
% % xmax = [3.56200500000000,33.0631965000000,35.6877540000000,9.06129750000000,1.26274350000000,66.8020020000000,181.665840000000,8.16812250000000,3.35976900000000,38.7809730000000,87.0764610000000,11.3072955000000,9062.55183300000,14.8805715000000,0.742906500000000];
% 
% distrpar=cell(M,1); for i=1:M; distrpar{i}=[xmin(i) xmax(i)]; end
% 
% % Load model output:
% % Set up parallel pool
% myCluster = parcluster('local');
% myCluster.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE')) / str2double(getenv('SLURM_CPUS_PER_TASK'));
% myCluster.JobStorageLocation = getenv('TMPDIR');
% myPool = parpool(myCluster, myCluster.NumWorkers);
% 
% 
% % Define name of simulink model
% mdlName     = 'ICN_with_BR_input_model4test_clusterParams_ICNtune_v8'%'ICN_with_BR_input_model4test_clusterParams_ICNnoLCN_v2';%'ICN_with_BR_input_model4test_clusterParams_HomogeneousICN_v2';
% 
% %% Step 3: Apply PAWN
% % (edit) these values
% N = 100; n = 6;
% DistrFun = 'unif' ;
% 
% 
% X = AAT_sampling('lhs',M,DistrFun,distrpar,N);
% 
% 
% % Feed X into model to get model output of interest Yu (CO, HR, MAP, etc)
% % Create an array of simulation input objects and specify the sweep value for each simulation
% simIn(1:length(X)) = Simulink.SimulationInput(mdlName);
% for idx = 1:length(X)
%     simIn(1,idx) = simIn(1,idx).setModelParameter('SaveTime', 'on', ...
%         'SaveOutput', 'on', ...
%         'TimeOut', 240);
% 
%     % Base model
%     simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(X(idx,1)), ...        % Set ICN parameters
%         [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(X(idx,2)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(X(idx,3)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(X(idx,4)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(X(idx,5)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(X(idx,6)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(X(idx,7)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(X(idx,8)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(X(idx,9)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(X(idx,10)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(X(idx,11)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(X(idx,12)), ...
%         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(X(idx,13)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(X(idx,14)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_BRgain'], 'Value', num2str(X(idx,15)), ...
%         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(X(idx,16)));
% 
%     % Homog ICN
% %         simIn(idx) = simIn(idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/ICN-TF/minval_ICN'], 'Value', num2str(X(idx,1)), ... % NActr
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/fmax_ICN'], 'Value', num2str(X(idx,2)), ... %deltaNActr), ...
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/midptICN'], 'Value', num2str(X(idx,3)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/kICN'], 'Value', num2str(X(idx,4)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/LCN_fevEmaxgain'], 'Value', num2str(X(idx,5)), ... % BR cell grp
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/LCN_BRgain'], 'Value', num2str(X(idx,6)), ... %(deltaBR), ...
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/LCN_feshgain'], 'Value', num2str(X(idx,7)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/ICN-TF/LCN_fevHRgain'], 'Value', num2str(X(idx,8)));
% 
% % No LCN
% % simIn(1,idx) = simIn(1,idx).setBlockParameter([mdlName '/Autonomic Nervous System/ICN/PN-NA/minval_NA_PN'], 'Value', num2str(X(idx,1)), ...        % Set ICN parameters
% %         [mdlName '/Autonomic Nervous System/ICN/PN-NA/fmax_NA_PN'], 'Value', num2str(X(idx,2)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-NA/midptNA_PN'], 'Value', num2str(X(idx,3)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-NA/kNA'], 'Value', num2str(X(idx,4)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/LCN/minval_LCN'], 'Value', num2str(X(idx,5)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/LCN/fmax_LCN'], 'Value', num2str(X(idx,6)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/LCN/midptLCN'], 'Value', num2str(X(idx,7)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/LCN/kLCN'], 'Value', num2str(X(idx,8)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/minval_DMV_PN'], 'Value', num2str(X(idx,9)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/fmax_DMV_PN'], 'Value', num2str(X(idx,10)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/midptDMV_PN'], 'Value', num2str(X(idx,11)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/kDMV_PN'], 'Value', num2str(X(idx,12)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/PN-DMV/tau_DMV_PN'], 'Value', num2str(X(idx,13)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_fevEmaxgain'], 'Value', num2str(X(idx,14)), ...
% %         [mdlName '/Autonomic Nervous System/ICN/LCN/LCN_feshgain'], 'Value', num2str(X(idx,15)));
%    
% end
% 
% % Run simulations
% simout = parsim(simIn);
% %simout = sim(simIn(1,2));
% 
% % Analysis for HR
% meanHR = zeros(length(X),1);
% for i = 1:length(X)
%     % Determine indices to use for 15 seconds at steady state
%     simTime         = simout(1,i).time;
%     tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
%     tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
%     %meanTime = mean(simout(1,i).time(tshortIdx));
%     
%     % Calculate avg HR over 15 seconds
%     meanHR(i,1)         = mean(simout(1,i).HR(tshortIdx));
% end
% 
% 
% Y = meanHR; % vector (2000,1)
% 
% % % Check Y for NaN and count number of NaN simulations
% % z = isnan(Y);
% % k = find(z); % indices of NaN values of Y
% % 
% % % Number of Nan values reomoved
% % numRemovedVals = length(k);
% % sprintf('Number of NaN simulation is %1.0f',numRemovedVals)
% % 
% % % Remove NaN values from X and Y
% % for i = length(k):-1:1
% %     Y(k(i)) = [];
% %     X(k(i),:) = [];
% % end
% 
% size(X)
% size(Y)
% 
% 
% [KS,KS_dummy,YF,Fu,Fc,YY2,xc2,NC,XX2,YU,idx_bootstrap] = pawn_ks_givendata(X,Y,n);
% 
% %% Plotting cdf and KS
% col=gray(n+1); % Color for colormap, (n+1) so that the last line is not white
% map = colormap(col(1:n,:));% Remove white color
% i =1; % input factor under investigation
% fs = 22; ms = 10 ; % font size and marker size
% 
% % (g): plot CDFs
% subplot(1,2,1)
% for k=1:n
%     plot(YF,Fc{i,k},'Color',map(k,:),'LineWidth',3.5); hold on
% end
% plot(YF,Fu,'r','LineWidth',3.5)
% set(gcf, 'Position',  [100, 100, 1200, 500])
% set(gca,'FontSize',fs); xlabel('y'); ylabel('cdf'); 
% % (h): plot KS
% subplot(1,2,2)
% for k=1:n
%     plot(xc2{i}(k),KS(k,i),'ok','MarkerFaceColor',map(k,:),'MarkerSize',ms); hold on
% end
% %set(gca,'FontSize',fs,'XTick',round(xc2{i},1),'YLim',[0,1],'XGrid','on')
% xlabel('x_1'); ylabel('KS'); 
% 
% saveas(gcf,'cdfandKS_HR_Homog.png')
% 
% %% Barchart of KS stats
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compute PAWN sensitivity indices for different choices of N
% % and different aggregation statistics
% % (this will take some computing time!!!)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% n = 6 ; % reduced from 10
% %N = N - numRemovedVals;
% N_red = [ 2 40 60 80 ] ; % size of reduced sample for approximation of PAWN indices
% 
% nboot = 5  ; %reduced from 50
% alfa  = 0.05 ;
% 
% PAWN_median   = nan(length(N_red),M); PAWN_median_lb = PAWN_median; PAWN_median_ub = PAWN_median;
% PAWN_mean     = nan(length(N_red),M); PAWN_mean_lb = PAWN_mean; PAWN_mean_ub = PAWN_mean;
% PAWN_max      = nan(length(N_red),M); PAWN_max_lb = PAWN_max; PAWN_max_ub = PAWN_max;
% KS_dummy_mean = nan(1,length(N_red));
% 
% for r=1:length(N_red) % WARNING: running this for loop may take some time!
%     r
%     idx = randperm(N,N_red(r)) ;
%     [KS_median,KS_mean,KS_max,KS_dummy] = pawn_indices_givendata(X(idx,:),Y(idx),n,nboot) ;
%     
%     % Take statistics across bootstrap resamples:
%     PAWN_median(r,:) = mean(KS_median) ;
%     median_lb = sort(KS_median) ; median_lb = median_lb(max(1,round(nboot*alfa/2)),:);
%     median_ub = sort(KS_median) ; median_ub = median_ub(round(nboot*(1-alfa/2)),:);
%     PAWN_median_lb(r,:) = median_lb ;
%     PAWN_median_ub(r,:) = median_ub ;
%     %
%     PAWN_mean(r,:)   = mean(KS_mean)  ;
%     mean_lb = sort(KS_mean) ; mean_lb = mean_lb(max(1,round(nboot*alfa/2)),:);
%     mean_ub = sort(KS_mean) ; mean_ub = mean_ub(round(nboot*(1-alfa/2)),:);
%     PAWN_mean_lb(r,:) = mean_lb ;
%     PAWN_mean_ub(r,:) = mean_ub ;
%     %
%     PAWN_max(r,:)   = mean(KS_max)  ;
%     max_lb = sort(KS_max) ; max_lb = max_lb(max(1,round(nboot*alfa/2)),:);
%     max_ub = sort(KS_max) ; max_ub = max_ub(round(nboot*(1-alfa/2)),:);
%     PAWN_max_lb(r,:) = max_lb ;
%     PAWN_max_ub(r,:) = max_ub ;
%     %
%     KS_dummy_mean(r) = mean(KS_dummy)  ;
%     %
%     fprintf('iteration %d (N_red=%d) of %d completed \n',r,N_red(r),length(N_red))
% 
% end
% % 
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Figure 4: Bar plot of PAWN indices for given (N,n) [using stat=median] 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% [tmp,ranking] = sort(-PAWN_median(end,:));
% %X_Labels_swat_ranked=X_Labels_swat(ranking);
% 
% 
% hfig= figure; 
% fs = 24 ;
% b = bar(PAWN_median(end,ranking)) ;
% set(b,'FaceColor',[126 126 126]/256)
% axis([0,M+1,0,0.6])
% hold on
% plot([1:M],ones(1,M)*KS_dummy_mean(end),'r','LineWidth',2) % plots red line
% %text(3,0.4,['N = ' int2str(N_red(end)) ', n = ' num2str(n) ],'FontSize',fs)
% set(gca,'FontSize',fs)
% for i=1:M; Xlab{i}=num2str(ranking(i)); end
% set(gca,'XTick',[1:M],'XTickLabel',Xlab,'YTick',[0,0.2,0.4],'FontSize',16)
% %for i=1:M; text(i-0.5,PAWN_median_ub(end,ranking(i))+0.02,num2str(ranking(i)),'FontSize',20); end
% ylabel('median(KS)','FontSize',fs)
% xlabel('input factors','FontSize',fs)
% for i = 1:M % connect upper and lower bound with a line
%     plot([i i],[PAWN_median_lb(end,ranking(i)) PAWN_median_ub(end,ranking(i))],'k','LineWidth',1.5)
% end
% set(hfig, 'Position', [0 0 1600 400])
% saveas(gcf,'newPAWNbar_HR_Homog.png')

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
saveas(gcf,'newPAWNbarCompare_HR_916030.png')

save('PAWN_MAP_noLCN.mat','PAWN_median','KS_dummy_mean','PAWN_median_lb','PAWN_median_ub')



%%
% % Estimate unconditional and conditional CDFs:
% [ YF, Fu, Fc  ] = pawn_cdfs(Yu,YY) ; % YY = (M,n) cell; each cell (NC,1); Yu = (NU,1)
% 
% % Plot CDFs:
% figure
% for i=1:M
%     subplot(1,M,i)
%     pawn_plot_cdf(YF, Fu, Fc(i,:),[],'y (HR)')
% end
% 
% % Save plot
% saveas(gcf,'cdfplot.png')
% 
% % Further analyze CDF of one input:
% i = 3 ;
% figure;
% pawn_plot_cdf(YF, Fu, Fc(i,:),xc{i},'y (HR)',labelparams{i}) % same
% % function as before but exploiting more optional input arguments
% 
% % Save plot
% saveas(gcf,'cdfplotNAPNmidpt.png')
% 
% % Compute KS statistics:
% KS = pawn_ks(YF,Fu,Fc) ;
% 
% % Plot KS statistics:
% figure
% for i=1:M
%     subplot(1,M,i)
%     pawn_plot_kstest(KS(:,i),NC,NU,0.05,xc{i},labelparams{i})
% end
% 
% % Save plot
% saveas(gcf,'KSstat.png')
% 
% % Compute PAWN index by taking a statistic of KSs (e.g. max):
% Pi = max(KS);
% 
% % Plot:
% figure
% boxplot1(Pi,labelparams)
% 
% % Save plot
% saveas(gcf,'boxplotPAWNidx.png')
% 
% % Use bootstrapping to assess robustness of PAWN indices:
% stat = 'max' ; % statistic to be applied to KSs
% Nboot = 100  ; % number of boostrap resamples
% [ T_m, T_lb, T_ub ] = pawn_indices(Yu,YY,stat,[],Nboot);
% 
% % Plot:
% figure; boxplot1(T_m,labelparams,[],T_lb,T_ub)
% 
% % Save plot
% saveas(gcf,'boxplotMax.png')
% 
% % Convergence analysis:
% stat = 'max' ; % statistic to be applied to KSs
% NCb = [ NC/10 NC/2 NC ] ;
% NUb = [ NU/10 NU/2 NU ] ;
% 
% [ T_m_n, T_lb_n, T_ub_n ] = pawn_convergence( Yu, YY, stat, NUb, NCb,[],Nboot );
% NN = NUb+n*NCb ;
% figure; plot_convergence(T_m_n,NN,T_lb_n,T_ub_n,[],'no of evals',[],labelparams)
% 
% % Save plot
% saveas(gcf,'convergenceplot.png')

%% Step 4: Apply PAWN to sub-region of the output range

% % Compute the PAWN index over a sub-range of the output distribution, for
% % instance only output values above a given threshold
% thres = 50 ;
% [ T_m2, T_lb2, T_ub2 ]= pawn_indices( Yu, YY, stat,[], Nboot,[],'above',thres ) ;
% 
% % Plot:
% figure; boxplot1(T_m2,labelparams,[],T_lb2,T_ub2)
% 
% % Save plot
% saveas(gcf,'boxplotSubregion.png')

%% Exit code
close_system
delete(myPool);
exit