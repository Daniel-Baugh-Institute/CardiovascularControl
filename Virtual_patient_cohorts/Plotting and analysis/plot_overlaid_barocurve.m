function plot_overlaid_barocurve(postIRstore_filename,color,postIRbarocurve_filename)
% Plot overlaid baroreflex curve with Osculati 1990 MI data
addpath 'C:\Users\mmgee\Box\Michelle-Gee\Research\MI model'
addpath(genpath('C:\Users\mmgee\AppData\Local\Temp\Mxt231\RemoteFiles'))
addpath(genpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/ACC2024-02-02-2024/'))
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/UnitTests/')
addpath('/lustre/ogunnaike/users/2420/matlab_example/matlab_slurm/Mastitskaya2012/')


%% Fit Osculati 1990 mean data
% load data
neckChamberPressure = [8	13.2	18.5	23.7	30	37.5]; % mm hg
RRbase = 777.8; % +/- 34.7 msec
deltaRR = [46.50360279	64.79801636	89.04640376	126.5452133	144.0065953	171.9898378];
SEM = [7.550265772	10.06205045	13.11093148	24.29514847	24.91603345	24.7516633]; % msec
num = 12;
SD = SEM.* sqrt(num); % standard deviation
CI95 = 1.96.* SD;
RR_Osculati = RRbase + deltaRR;

%% Plot
% data
x = neckChamberPressure;
y = RR_Osculati;
errors = CI95;


%% Plot accepted baroreflex curves based on meeting baseline HR, MAP criteria and baroreflex curve criteria (1, 2, 6)
% loop to load mat files for ['postIR_cardiac_' num2str(i) '_barocurve.png.mat']-- gives BPvalues, HR-- plot in a loop over Osculati 1990 data
ms = 8;
numPatients = 59;
numSampleSets = 100;
osculatiIdx = 4:10;
load(postIRstore_filename)
count = 0;

figure;
hold on;

% postIRstore_filename = 'postIR_baroreceptor_06182024';
for i = 1:numPatients
    filename = [postIRbarocurve_filename '_' num2str(i) '.mat'];
    load(filename,'HR','BPvalues')
    for j = 1:numSampleSets
        if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
            RR = 1000.*60./HR(j,:);
            if count > 0
                plot(BPvalues,RR,'-','Color',color,'HandleVisibility','off')
            else
                plot(BPvalues,RR,'-','Color',color,'HandleVisibility','on')
                count = count + 1;
            end
        end
            
    end
end

% % for NADMV/NTS overlay
% postIRstore_filename = 'postIR_NTS_06182024';
% load(postIRstore_filename)
% green = [0.4660 0.6740 0.1880];
% color = green;

% count = 0;
% for i = 1:numPatients
%     filename = [postIRstore_filename '_' num2str(i) '.mat'];
%     load(filename,'HR','BPvalues')
%     for j = 1:numSampleSets
%         if postIRstore(i).crit(j,1) == 1 && postIRstore(i).crit(j,2) == 1 && postIRstore(i).crit(j,6) == 1
%             RR = 1000.*60./HR(j,:);
%             if count > 0
%                 plot(BPvalues,RR,'-','Color',color,'HandleVisibility','off')
%             else
%                 plot(BPvalues,RR,'-','Color',color,'HandleVisibility','on')
%                 count = count + 1;
%             end
%         end
% 
%     end
% end

errorbar(x, y, errors,'ro', 'MarkerFaceColor','r','MarkerSize',ms);

xlabel('Neck chamber pressure (mm Hg)');
ylabel('RR interval (ms)');
% ttl = title('E');
% ttl.Units = 'Normalize';
% ttl.Position(1) = -0.33; % use negative values (ie, -0.1) to move further left
% ttl.HorizontalAlignment = 'left';
legend('Model', 'Osculati 1990') %'NTS model',
set(gca,'FontSize',14)

set(gcf,'Position',[10,10,400,400])
hold off;
saveas(gcf,[postIRstore_filename '.png'])
saveas(gcf,[postIRstore_filename '.fig'])
end