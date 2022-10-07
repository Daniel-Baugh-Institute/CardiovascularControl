function CO_SV = fx_CO_SV(tRange, time, For, Fol, Tperiod)
% Adapted from Park 2020
%% fx calculates cardiac output (CO) & stroke volume for each pulse of the heart and returns the CO for [Rventr, Lventr] and SV for [Rventr, Lventr]
%% return structure is a 4-column matrix [Right-CO, Left-CO, Right-SV, Left-SV]

index           = find(time >= tRange(1) & time <= tRange(2));
pkThreshold     = 25;  % lower limit (mL) threshold that determines if potential peak is representative of heart beat

[For_pks, Rpk_index] = findpeaks(For(index), 'MinPeakHeight', pkThreshold,  'MinPeakProminence', 5);   % returns peak values and indices of peaks
[Fol_pks, Lpk_index] = findpeaks(Fol(index), 'MinPeakHeight', pkThreshold,  'MinPeakProminence', 5);

% % remove any residual shoulders mistaken for peaks
% threshold               = 25;   % vol limit to be considered a peak    
% For_idx2rm              = find(For_pks < threshold);
% For_pks(For_idx2rm)     = []; 
% Rpk_index(For_idx2rm)   = [];
% 
% Fol_idx2rm              = find(Fol_pks < threshold);
% Fol_pks(Fol_idx2rm)     = [];
% Lpk_index(Fol_idx2rm)   = [];


% Right ventricle
Rpk_times = time(index);            % get time values 
Rpk_times = Rpk_times(Rpk_index);   % identify time at which peaks occur 
Rpk_times = Rpk_times(2:(end-1));   % exclude first and last peaks to avoid issues with mismatch dimensions later on

% Left ventricle
Lpk_times = time(index);
Lpk_times = Lpk_times(Lpk_index);
Lpk_times = Lpk_times(2:(end-1));  

% find index of times for each beat within tRange
RH_per          = Tperiod(index);
RH_per          = RH_per(Rpk_index); % Tperiod at the peak of For
RH_per          = RH_per(2:(end-1));    % exclude first and last time points
RHpkEndTimes    = Rpk_times + RH_per/2; % end of stroke is calculated as time of peak outflow from ventricle plus half of heart period

[Rres1, Rres2]  = histcounts(RHpkEndTimes, time(index)); % determine index of actual times that are closest to end peak times
% [Rres1, Rres2]  = hist(RHpkEndTimes, time(index)); 

cumRV           = cumtrapz(time(index), For(index));
RSV             = cumRV(find(Rres1 >= 1)+3); % right ventricular stroke volume

if length(RSV)  == 1
    RCO     = RSV/RH_per;
else
    RSV     = diff(RSV);
    if length(RSV) == length(RH_per(2:end-1)) % edited from script provided by Park 2020, which did not have this if statement and skipped directly to else statement, which would produce an NaN if the vectors were not the same length
        RCO     = RSV./RH_per(2:end-1);
        disp('Warning: RCO length changed')
    elseif length(RSV) == length(RH_per(2:end-2))
        RCO     = RSV./RH_per(2:end-2);
        disp('Warning: RCO length changed')
    elseif length(RSV) == length(RH_per(2:end))
        RCO     = RSV./RH_per(2:end);
    else
        disp('Warning: Right Cardiac Output issue')
        RCO = NaN;
    end
end

LH_per          = Tperiod(index);
LH_per          = LH_per(Lpk_index);
LH_per          = LH_per(2:(end-1));    % exclude first and last time points
LHpkEndTimes    = Lpk_times + LH_per/2;
% [Lres1, Lres2]  = histcounts(LHpkEndTimes, time(index)); % determine index of actual times that are closest to end peak times
[Lres1, Lres2]  = hist(LHpkEndTimes, time(index)); % determine index of actual times that are closest to end peak times


cumLV           = cumtrapz(time(index), Fol(index));
LSV             = cumLV(find(Lres1 >= 1)+3); 
% LSV(2:end)      = diff(LSV);

if length(LSV)  == 1
    LCO     = LSV/LH_per;
else
    LSV     = diff(LSV);
    if length(LSV) == length(LH_per(2:end-1)) % edited from script provided by Park 2020, which did not have this if statement and skipped directly to else statement, which would produce an NaN if the vectors were not the same length
        % Mismatching or matching indices doesn't seem to change output of
        % calculation
        LCO     = LSV./LH_per(2:end-1);
        disp('Warning: LCO length changed')
    elseif length(LSV) == length(LH_per(2:end-2))
        LCO     = LSV./LH_per(2:end-2);
        disp('Warning: LCO length changed')
    elseif length(LSV) == length(LH_per(2:end))
        LCO     = LSV./LH_per(2:end);
    else
        disp('Warning: Right Cardiac Output issue')
        LCO     = NaN;
    end
end
CO_SV = struct('RCO', RCO, 'RSV', RSV, 'LCO', LCO,  'LSV', LSV);

