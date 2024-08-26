function sigmoid_compare(part,sig_values,filename,legend_entries)
% Inputs:
% part: NTS = 1, NA/DMV = 2, ICN = 3, baroreflex = 4
% sig_values: 1xn vector of new sigmoid values
% filename: string for filename, do not include suffix
% legend_entries: cell array with strings for legend entries

% parameters
% ICN
ICNparams = [0.857242 12.624676 13.811104 3.230162; 1.988567 18.375719 521.217197 3.003231; 2.415242 17.718990 14.183806 13.356069];% PNDMVdelay 3.329861 2.661685 5.642977 0.066794];
% NA/DMV
NAparams = [4.88, 15.78, 59.83, 23; 0.61, 11, 12.81, 7; 2.5901, 6.66, 42.91, 33.5]; %  NA, NActr, DMV
% NTS
NTSparams = [0.30, 21.50, 37.07, 21; 0.45, 28.33, 10.2, 7; 2.75, 31.57, 11.13, 2]; %  BR, CPR, LSR; fmin, fmax, fmid, k
% baroreflex
ka = 11.758;
P_n = 92;
f_min = 2.52;
f_max = 47.78;



if part == 1
    params = NTSparams;
elseif part == 2
    params = NAparams;
elseif part == 3
    params = ICNparams;
else
    params = ka;
end

% pre-ischemia sigmoids
x = 0:0.1:1000;
if part < 4
    for i = 1:3
        fmin = params(i,1);
        fmax = params(i,2);
        fmid = params(i,3);
        k = params(i,4);
        y_PI(i,:) = (fmin-fmax)./(1+(x./fmid).^k) + fmax;
    end
else
    y_PI = (f_min + f_max.*exp((x-P_n)./ka))./(1+exp((x-P_n)./ka));
end

% end-ischemia sigmoids
if part < 4
    for i = 1:3
        fmin = sig_values(i,1);
        fmax = sig_values(i,2);
        fmid = sig_values(i,3);
        k = sig_values(i,4);
        y_EI(i,:) = (fmin-fmax)./(1+(x./fmid).^k) + fmax;
    end
else
    y_EI = (f_min + f_max.*exp((x-P_n)./sig_values))./(1+exp((x-P_n)./sig_values));
end

% plot formatting
% plot bounds
xlim_NTS = [10,50;0,20;0,50];
xlim_NA = [40,80; 0,25; 30,48];
xlim_ICN = [0,25;0,1000; 5,30];
lims = [xlim_NTS; xlim_NA; xlim_ICN];
xlim_ka = [40 160];


[num_sig,num_params] = size(sig_values);
rows = 1;
cols = num_sig;
fs = 16;
lw = 4;
markerColors = {'b','g','m','r','c'};

if part < 4
    for i = 1:3
        subplot(rows,cols,i)
        plot(x,y_PI(i,:),'LineWidth',lw,'Color',markerColors{1},'DisplayName',legend_entries{1});
        hold on
        plot(x,y_EI(i,:),'LineWidth',lw,'Color',markerColors{2},'DisplayName',legend_entries{2});
        xlabel('Net input activity (Hz)')
        ylabel('Output firing rate (Hz)')
        set(gca,'FontSize',fs)
        hold off
        if i == 3
            legend('Location','Northwest')%legend_entries,
        end
        if part == 1
            limits = lims(i,:);
        elseif part == 2
            limits = lims(i+3,:);
        else
            limits = lims(i+6,:);
        end
        xlim(limits)
    end
    set(gcf, 'Position',  [10, 10, 1800, 400])
else
    size(y_PI)
    size(y_EI)
    plot(x,y_PI,'LineWidth',lw,'Color',markerColors{1},'DisplayName',legend_entries{1})
    hold on
    plot(x,y_EI,'LineWidth',lw,'Color',markerColors{2},'DisplayName',legend_entries{2})
    xlabel('Pressure (mm Hg)')
    ylabel('Output firing rate (Hz)')
    xlim(xlim_ka)
    legend('Location','Northwest')
    set(gcf, 'Position',  [10, 10, 600, 400])
    set(gca,'FontSize',fs)
end




file_suffix = '.png';
png_filename = [filename file_suffix];
saveas(gcf,png_filename)
end

