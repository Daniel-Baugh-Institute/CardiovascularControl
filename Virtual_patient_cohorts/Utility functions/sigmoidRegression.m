function [A1,A2,A3,A4] = sigmoidRegression(ECSP,RR)
% fit baroreflex curve data to a sigmoidal curve based on Seredynski 2021,
% Fadel 2003 function
% Response variable = A1{1+exp[A2(ECSP-A3)]}^-1 +A4

% Inputs: 
    % ECSP: estimated carotid sinus pressure (neck chamber pressure + MAP)
    % RR: RR interval
% Outputs: fitted regression parameters, plot of data and fitted curve
    % A1: response range (max - min)
    % A2: gain coefficient
    % A3: centering point
    % A4: minimum response
    
% References:
% Fadel: Recent insights into carotid baroreflex function in humans using
% the variable pressure neck chamber
% Seredynski: Neck chamber technique revisited: low-noise device delivering
% negative and positive pressure and enabling concomitant carotid artery
% imaging with ultrasonography


% Define Start points, fit-function and fit curve
x0 = [0.2 0.01 100 0.6]; % MAP initial guess
x = ECSP;
fitfun = fittype( @(A1,A2,A3,A4,x) A1.*(1+exp(A2.*(x-A3))).^(-1) + A4 );
% Set the fit options
% options = fitoptions('Method', 'NonlinearLeastSquares', ...
%                      'Lower', [0 -10 0 0], ...
%                      'Upper', [100 0 200 5], ...
%                      'StartPoint', x0); % bounding makes the fit
%                      signficantly worse

options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', x0);%     
%                      

[fitted_curve,gof] = fit(ECSP,RR,fitfun,options);

% Save the coeffiecient values for a,b,c and d in a vector
coeffvals = coeffvalues(fitted_curve);
A1 = coeffvals(1);
A2 = coeffvals(2);
A3 = coeffvals(3);
A4 = coeffvals(4);

% Plot results
fs = 16;
lw = 2;


ECSP_plot = sort(ECSP); % sort by ascending order for plotting fitted curve
% scatter(ECSP, RR, 'r+')
hold on
RR_plot = fitted_curve(ECSP_plot);
plot(ECSP_plot,RR_plot,'b-','LineWidth',2)
% ylim([60 120])
xlabel('Estimated carotid sinus pressure (mm Hg)')
ylabel('RR interval (s)')
set(gca, 'FontSize',fs)
% legend('Data for 5 individuals','Fit','Location','NorthWest')
% hold off
% saveas(gcf,'fittedEnsembleSeredynski.png')
% save('sigmoidRegression.mat','ECSP_plot','RR_plot')
end