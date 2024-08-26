function [SDRR, RMSSD, pNN50] = HRVmetrics_timeDomain(simOut,trange)
% Input:
%   simOut: 1xn vector of simulation outputs
% Calculates
%     standard deviation of RR interval (SDRR)

%     root mean square of successive differences (RMSSD): "sensitive to
%     high-frequency heart period fluctuations in the respiratory frequency
%     range and has been used as an index of vagal cardiac control"
%     (Berntson et al 2005)

%     Low frequency to high frequency ratio

numMdls = length(simOut);
SDRR = zeros(numMdls,1);
RMSSD = zeros(numMdls,1);
pRMS = zeros(numMdls,1);
powbp_total = zeros(numMdls,1);
powbp_vlf = zeros(numMdls,1);
powbp_lf = zeros(numMdls,1);
powbp_hf = zeros(numMdls,1);
LFHFratio = zeros(numMdls,1);

for i = 1:numMdls
    % Calculate RR intervals
    [RRinterval, RRtimes] = RRinterval_calc(simOut,trange);
    RRinterval_msec = RRinterval * 1000;

    % SDRR
    SDRR(i,1) = std(RRinterval);

    % RMSSD
    RMSSD(i,1) = sqrt(mean((diff(RRinterval)).^2));

    % pNN50 
    RRdiff = diff(RRinterval_msec);
    NN50 = length(find(RRdiff > 50));
    pNN50 = NN50 / length(RRdiff);

end

end