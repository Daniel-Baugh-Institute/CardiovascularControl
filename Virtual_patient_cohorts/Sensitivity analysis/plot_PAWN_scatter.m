function plot_PAWN_scatter(MAPfile,HRfile,outFilename)
% Scatter plot of MAP vs HR sensitivity index
% Inputs: MAPfile - File containing MAP sensitivity data
    %         HRfile - File containing HR sensitivity data
    %         outFilename - Output filename for the plot
% Output: Plot with error bars

% load data
load(MAPfile)
PAWN_median_MAP = PAWN_median(4,:);
PAWN_median_ub_MAP = PAWN_median_ub(4,:) - PAWN_median(4,:);
PAWN_median_lb_MAP = PAWN_median(4,:) - PAWN_median_lb(4,:);
KS_dummy_mean_MAP = KS_dummy_mean(4);
load(HRfile)
PAWN_median_HR = PAWN_median(4,:);
PAWN_median_ub_HR = PAWN_median_ub(4,:) - PAWN_median(4,:);
PAWN_median_lb_HR = PAWN_median(4,:) - PAWN_median_lb(4,:);
KS_dummy_mean_HR = KS_dummy_mean(4);

% plot with error bars
figure;
    hold on;
    errorbar(PAWN_median_MAP, PAWN_median_HR, PAWN_median_lb_MAP, PAWN_median_ub_MAP, PAWN_median_lb_HR, PAWN_median_ub_HR, ...
        'o','Color',[0.7 0.7 0.7],'MarkerFaceColor','b','LineWidth',2,...
        'MarkerEdgeColor','none','MarkerSize',8);
    plot([0 0.3],[KS_dummy_mean_HR KS_dummy_mean_HR],'r--')
    plot([KS_dummy_mean_MAP KS_dummy_mean_MAP],[0 0.25],'r--')
    xlabel('Mean arterial pressure sensitivity index');
    ylabel('Heart rate sensitivity index');
    set(gca,'FontSize',16)
    set(gcf,'Position',[10 10 600 500])
    hold off;

    % Save plot to file
    if ~isempty(outFilename)
        saveas(gcf, outFilename);
    end

    % Print the indices of variables with high sensitivity index
    idxToPrint = find(PAWN_median_MAP >= 0.078);
    disp(idxToPrint)
end