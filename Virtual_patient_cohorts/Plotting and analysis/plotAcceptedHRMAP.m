function plotAcceptedHRMAP(preIRHRMAP,postIRHRMAP,numMdlsByType,saveFilename)
% Phase plot of HR and MAP for accepted pre and post IR models to see
% if there's any clustering\
% Inputs:
%   PreIRHRMAP: first row HR, 2nd row MAP
%   postIRHRMAP
%   numMdlsByType: vector with number of preIR, postIR cardiac, postIR BR, postIR
%   NTS, postIR NADMV, postIR ICN, postIR all models
%   saveFilename: Filename to save the figure

% generate indices for plotting
idx2 = cumsum(numMdlsByType(2:end));
idx1 = idx2-numMdlsByType(2:end)+1;



% plot formatting
purple = [0.4940 0.1840 0.5560];
goldenrod = [0.9290 0.6940 0.1250];
green = [0.4660 0.6740 0.1880];
magenta = [ 1 0 1];
cyan = [0 1 1];
brown = [ 0.7686    0.6000    0.4863];
blue = [0 0 1];
grey = [0.9020    0.9020    0.9020];
greyOutline = [0.6510    0.6510    0.6510];
black = [0 0 0];
white = [1 1 1];
greyblue =[0.8510    0.9647    0.9882];
color = [purple; goldenrod; green; magenta; cyan; brown];


figure;
plot(preIRHRMAP(1,:),preIRHRMAP(2,:),'p','MarkerEdgeColor',greyOutline,'MarkerFaceColor','b');
hold on
for i = 1:length(numMdlsByType)-1
    plot(postIRHRMAP(1,idx1(i):idx2(i)),postIRHRMAP(2,idx1(i):idx2(i)),'o','Color',color(i,:),'MarkerFaceColor',color(i,:));
end
legend('Pre-IR','Cardiac','Baroreceptors','NTS','NA+DMV','ICN','All','Location','northeastoutside')
ylabel('MAP (mm Hg)');
xlabel('HR (bpm)')
set(gca,'FontSize',15)
set(gcf,'Position',[100, 100, 700,400])

% Save the figure
if ~isempty(saveFilename)
    saveas(gcf, saveFilename);
end
end

