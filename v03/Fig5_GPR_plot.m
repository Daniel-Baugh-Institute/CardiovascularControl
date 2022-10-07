% Closed-loop modeling of intrinsic cardiac nervous system contributions to
% respiratory sinus arrhythmia 
% Michelle Gee 
% October 3, 2022

% Script to produce Figure 5. Script to run simulations can be found in
% Fig5_GPR_cluster.m, which runs on a high performance computing cluster
% (University of Delaware's DARWIN computing cluster was used in this
% case).

% Load data
close all; clear;
load('RSA_GPR_kRSA_09_30_22.mat')

% define kRSA values
LowerBound = 0;
UpperBound = 1;
kRSA = LowerBound:0.1:UpperBound;

% plot formatting
fs = 18; % font size

figure(1)
plot(kRSA,-log(var_base),'bo','MarkerFaceColor','b')
xlabel('kRSA value')
ylabel('Likelihood')
xlim([0 0.7])

set(gca,'FontSize',fs)
saveas(gcf,'Fig5_kRSA_likelihood.png')