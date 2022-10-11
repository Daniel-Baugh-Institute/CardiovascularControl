%% Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
% Michelle Gee October 11, 2022

% PURPOSE Analyzing data from Rajendran et al. 2019 data on ICN neuron
% recordings in response to VNS in pig to extract average firing frequency
% before, during, and after VNS.


% REFERENCES Data set from Pennsieve:
% https://app.pennsieve.io/N:organization:618e8dd9-f8d2-4dc4-9abb-c6aaab2e78a0/datasets/N:dataset:0d9454ca-43d9-4fdb-ac79-01c3574e8565/files/N:collection:063f029d-ebf1-49fc-803c-1c618f314ba4

% Protocol: https://www.protocols.io/view/pig-icn-recording-2jugcnw?step=4
% Protocol also summarized in main text

% Rajendran, Pradeep S., Keijiro Nakamura, Olujimi A. Ajijola, Marmar
% Vaseghi, J. Andrew Armour, Jeffrey L. Ardell, and Kalyanam Shivkumar.
% 2016. “Myocardial Infarction Induces Structural and Functional
% Remodelling of the Intrinsic Cardiac Nervous System.” The Journal of
% Physiology 594 (2): 321–41.

% Rajendran, Pradeep, Marmar Vaseghi, and Jeffrey Ardell. 2019. “Functional
% Recordings from the Pig Intrinsic Cardiac Nervous System (ICN).” SPARC
% Consortium. https://doi.org/10.26275/OWRI-MPSX.

% Shin, Hyun-Chool, Vikram Aggarwal, Soumyadipta Acharya, Marc H. Schieber,
% and Nitish V. Thakor. 2010. “Neural Decoding of Finger Movements Using
% Skellam-Based Maximum-Likelihood Decoding.” IEEE Transactions on
% Biomedical Engineering 57 (3): 754–60.

% PSEUDO CODE
%  identify time intervals for baseline, stim, recovery create vectors for
%  each of those periods. Use the skellam distribution (Rajendran et al.
%  2016, Shin et al. 2010) to calculate firing frequency repeat for each
%  struct (each ICN neuron). This must be done separately because each
%  struct is a different length

clear; close all;
%% load data
load('Pig013_ICNS15_Matlab.mat')

%% Estimate average firing frequency using poisson distribution and determine significance of differences in lambda using skellam distribution
 
 % time intervals for baseline, stim, recovery (based on files from
 % Rajendran 2019)
 lcvBaseTime = [3052.2 3112.2];
 lcvStimTime = [3112.2 3171.1];
 lcvRecovTime = [3171.1 3231.1];
 rcvBaseTime = [3457.1 3517.1];
 rcvStimTime = [3517.1 3576.0];
 rcvRecovTime = [3576.0 3636.0];
 
% Preallocate
lcvBaseff = zeros(1,7);
lcvStimff = zeros(1,7);
lcvRecovff = zeros(1,7);
rcvBaseff = zeros(1,7);
rcvStimff = zeros(1,7);
rcvRecovff = zeros(1,7);
p_lcv_Base_Stim = zeros(1,7);
p_lcv_Base_Recov = zeros(1,7);
p_rcv_Base_Stim = zeros(1,7);
p_rcv_Base_Recov = zeros(1,7);

% Extract data
% After consulting with the aurthors of Rajendran et al. 2019, the provided
% .mat file is organized so that recordings for each neuron are in their
% own struct (st_36_01, etc). The field times contains the times that a
% spike occured at. So to calculate firing frequency, we count the number
% of elements in times over a given time period.
 for i = 1:5
     % extract time for each struct into vector 'time01, time02, etc'
     tempTime = eval(['st_36_0' num2str(i) '.times;']);
     eval(['time0' num2str(i) '=tempTime;'])
     
     % Indices for each interval
     tlcvBaseTimeIdx       = find(tempTime>= lcvBaseTime(1) & tempTime <= lcvBaseTime(2));
     tlcvStimTimeIdx       = find(tempTime>= lcvStimTime(1) & tempTime <= lcvStimTime(2));
     tlcvRecovTimeIdx       = find(tempTime>= lcvRecovTime(1) & tempTime <= lcvRecovTime(2));
     trcvBaseTimeIdx       = find(tempTime>= rcvBaseTime(1) & tempTime <= rcvBaseTime(2));
     trcvStimTimeIdx       = find(tempTime>= rcvStimTime(1) & tempTime <= rcvStimTime(2));
     trcvRecovTimeIdx       = find(tempTime>= rcvRecovTime(1) & tempTime <= rcvRecovTime(2));
     
     % Vectors for each period
      %tempTime = eval(['st_36_0' num2str(i) '.values(:,1);']);
      lcvBase = tempTime(tlcvBaseTimeIdx);
      lcvStim = tempTime(tlcvStimTimeIdx);
      lcvRecov = tempTime(tlcvRecovTimeIdx);
      rcvBase = tempTime(trcvBaseTimeIdx);
      rcvStim = tempTime(trcvStimTimeIdx);
      rcvRecov = tempTime(trcvRecovTimeIdx);
      
      % Normalize so each time epoch starts at zero for fitting to poisson
      % distribution
      lcvBaseNorm = lcvBase - lcvBaseTime(1);
      lcvStimNorm = lcvStim - lcvStimTime(1);
      lcvRecovNorm = lcvRecov - lcvRecovTime(1);
      rcvBaseNorm = rcvBase - rcvBaseTime(1);
      rcvStimNorm = rcvStim - rcvStimTime(1);
      rcvRecovNorm = rcvRecov - rcvRecovTime(1);
      
      
      % Estimate lambdas and calculate probability of differenes in lambda
      [p_lcv_Base_Stim(i),lcvBaseff(i),lcvStimff(i)] = skellamProb(lcvBaseNorm,lcvStimNorm,60);
      [p_lcv_Base_Recov(i),lcvBaseff(i),lcvRecovff(i)] = skellamProb(lcvBaseNorm,lcvRecovNorm,60);
      [p_rcv_Base_Stim(i),rcvBaseff(i),rcvStimff(i)] = skellamProb(rcvBaseNorm,rcvStimNorm,60);
      [p_rcv_Base_Recov(i),rcvBaseff(i),rcvRecovff(i)] = skellamProb(rcvBaseNorm,rcvRecovNorm,60);
      
 end
 
 % Repeat for structs 8 and 9 BUT RE-LABEL AS CELLS 6 AND 7
 for i = 8:9
     % extract time for each struct into vector 'time01, time02, etc'
     tempTime = eval(['st_36_0' num2str(i) '.times;']);
     eval(['time0' num2str(i) '=tempTime;'])
     
     % Indices for each interval
     tlcvBaseTimeIdx       = find(tempTime>= lcvBaseTime(1) & tempTime <= lcvBaseTime(2));
     tlcvStimTimeIdx       = find(tempTime>= lcvStimTime(1) & tempTime <= lcvStimTime(2));
     tlcvRecovTimeIdx       = find(tempTime>= lcvRecovTime(1) & tempTime <= lcvRecovTime(2));
     trcvBaseTimeIdx       = find(tempTime>= rcvBaseTime(1) & tempTime <= rcvBaseTime(2));
     trcvStimTimeIdx       = find(tempTime>= rcvStimTime(1) & tempTime <= rcvStimTime(2));
     trcvRecovTimeIdx       = find(tempTime>= rcvRecovTime(1) & tempTime <= rcvRecovTime(2));
     
     % Vectors for each period
      %tempTime = eval(['st_36_0' num2str(i) '.values(:,1);']);
      lcvBase = tempTime(tlcvBaseTimeIdx);
      lcvStim = tempTime(tlcvStimTimeIdx);
      lcvRecov = tempTime(tlcvRecovTimeIdx);
      rcvBase = tempTime(trcvBaseTimeIdx);
      rcvStim = tempTime(trcvStimTimeIdx);
      rcvRecov = tempTime(trcvRecovTimeIdx);
      
      % Normalize so each time epoch starts at zero for fitting to poisson
      % distribution
      lcvBaseNorm = lcvBase - lcvBaseTime(1);
      lcvStimNorm = lcvStim - lcvStimTime(1);
      lcvRecovNorm = lcvRecov - lcvRecovTime(1);
      rcvBaseNorm = rcvBase - rcvBaseTime(1);
      rcvStimNorm = rcvStim - rcvStimTime(1);
      rcvRecovNorm = rcvRecov - rcvRecovTime(1);
      
      % Estimate lambdas and calculate probability of differenes in lambda
      [p_lcv_Base_Stim(i-2),lcvBaseff(i-2),lcvStimff(i-2)] = skellamProb(lcvBaseNorm,lcvStimNorm,60);
      [p_lcv_Base_Recov(i-2),lcvBaseff(i-2),lcvRecovff(i-2)] = skellamProb(lcvBaseNorm,lcvRecovNorm,60);
      [p_rcv_Base_Stim(i-2),rcvBaseff(i-2),rcvStimff(i-2)] = skellamProb(rcvBaseNorm,rcvStimNorm,60);
      [p_rcv_Base_Recov(i-2),rcvBaseff(i-2),rcvRecovff(i-2)] = skellamProb(rcvBaseNorm,rcvRecovNorm,60);
      
 end  
%% plot firing frequencies
figure(22)
subplot(2,3,1)
bar(lcvBaseff)
ylabel({'Left cervical vagus stimulation'; 'average firing frequency'})
xlabel('Cell Number')
title('Baseline ff')
ylim([0 4])

subplot(2,3,2)
bar(lcvStimff)
ylabel('average firing frequency')
xlabel('Cell Number')
title('Stimulation ff')
ylim([0 4])

subplot(2,3,3)
bar(lcvRecovff)
ylabel('average firing frequency')
xlabel('Cell Number')
title('Recovery ff')
ylim([0 4])

subplot(2,3,4)
bar(rcvBaseff)
ylabel({'Right cervical vagus stimulation'; 'average firing frequency'})
xlabel('Cell Number')
ylim([0 4])

subplot(2,3,5)
bar(rcvStimff)
ylabel('average firing frequency')
xlabel('Cell Number')
ylim([0 4])

subplot(2,3,6)
bar(rcvRecovff)
ylabel('average firing frequency')
xlabel('Cell Number')
ylim([0 4])

disp('Average Firing Frequency')
disp('                          Baseline         Stimulation         Recovery')
fprintf('Left Vagus Stimulation \t\t%1.2f \t\t\t%1.2f \t\t\t\t%1.2f\n', mean(lcvBaseff),mean(lcvStimff),mean(lcvRecovff))
fprintf('Right Vagus Stimulation \t%1.2f \t\t\t%1.2f \t\t\t\t%1.2f\n\n\n', mean(rcvBaseff),mean(rcvStimff),mean(rcvRecovff))

disp('Probability of difference in firing from baseline for left VNS')
disp('                Stimulation           Recovery')
fprintf('Neuron 1  \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(1),p_lcv_Base_Recov(1))
fprintf('Neuron 2 \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(2),p_lcv_Base_Recov(2))
fprintf('Neuron 3 \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(3),p_lcv_Base_Recov(3))
fprintf('Neuron 4 \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(4),p_lcv_Base_Recov(4))
fprintf('Neuron 5 \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(5),p_lcv_Base_Recov(5))
fprintf('Neuron 6 \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(6),p_lcv_Base_Recov(6))
fprintf('Neuron 7 \t\t\t%1.2f \t\t\t\t%1.2f\n', p_lcv_Base_Stim(7),p_lcv_Base_Recov(7))
