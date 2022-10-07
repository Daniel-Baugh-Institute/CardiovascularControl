function simRes = PhysOutputs_Gen(simOut)
% Adapted from Park 2020
% takes simulation output and returns systolic blood pressure, diastolic
% blood pressure, heart rate, cardiac output, and ejection fraction
simRes = zeros(1,6);

if isnan(simOut.time)
    simRes = NaN;

else
    % start & end times for calculating steady state values
    tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
    tlong  = [120 180];     % for CO and SV determination

    simTime         = simOut.time;
    tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
    tlongIdx        = find(simTime>= tlong(1)  & simTime <= tlong(2));


    Emaxlv_vals     = simOut.Emaxlv;
    Tperiod_vals    = simOut.Tperiod;
    HR_vals         = simOut.HR_TS;


    Rsp_vals        = simOut.Rsp;
    Rep_vals        = simOut.Rep;
    Rmp_vals        = simOut.Rmp;


    Vlv_vals        = simOut.Vlv;


    Psa_vals        = simOut.Psa;

    Flow_lv_vals    = simOut.Flow_lv;
    Flow_rv_vals    = simOut.Flow_rv;

    Pla_vals        = simOut.Pla;
    Plv_vals        = simOut.Plv;

    fevHR_vals      = simOut.fevHR;
    fevEmax_vals    = simOut.fevEmax;

    phi_vals        = simOut.Phi;

    try
        CO_SV = fx_CO_SV(tlong, simTime, Flow_lv_vals, Flow_rv_vals, Tperiod_vals);
    catch
        warning('error in CO/SV calculation likely due to parameters used');
        CO_SV = struct('RCO', NaN, 'RSV', NaN, 'LCO', NaN,  'LSV', NaN);
    end

    EDVthreshold    = 10; % end diastolic volume min threshold
    try
        [Vlv_pks, Vlv_index]    = findpeaks(Vlv_vals(tshortIdx), 'MinPeakHeight', EDVthreshold,  'MinPeakProminence', 5);   % returns peak values and indices of peaks
    catch
        warning('error in left ventricular volume calculation')
        Vlv_pks     = 999;
    end

    meanVlvPks      = mean(Vlv_pks);
%     temp_EjFr       = mean(CO_SV.LSV)/mean(Vlv_pks);
    SysThreshold    = 10;

    try
        [Psys_pks, Psys_index]  = findpeaks(Psa_vals(tlongIdx),    'MinPeakHeight', SysThreshold, 'MinPeakProminence', 2); %mean peak prominence 5 
        [Pdia_pks, Pdia_index]  = findpeaks(-1*Psa_vals(tlongIdx), 'MinPeakProminence', 2);

        if length(Psys_pks) ~= length(Pdia_pks)
            minLength   = min([length(Psys_pks), length(Pdia_pks)]);
            Psys_pks    = Psys_pks(1:minLength);
            Pdia_pks    = Pdia_pks(1:minLength);
        end
    catch
        Psys_pks    = 999*ones(10, 1);  % arbitrarily chosen vector length of 10
        Pdia_pks    = 999*ones(10, 1);  % arbitrarily chosen vector length of 10
    end

    % identify appropriate indices for Left Ventricualr End Systolic Volume
    % (ESV) and End Systolic Pressure (ESP)
    Vlv_vals_tshort     = Vlv_vals(tshortIdx);
    phi_vals_tshort     = phi_vals(tshortIdx);
    Vlv_idx             = find(Vlv_vals_tshort <= (meanVlvPks-10));
    phi_idx             = find(phi_vals_tshort <= 0.2);
    int_idx             = intersect(Vlv_idx, phi_idx);
    ESV_idx             = int_idx(find(diff(int_idx)> 5) + 1);
    ESV_vals            = Vlv_vals_tshort(ESV_idx);

    Plv_vals_tshort     = Plv_vals(tshortIdx);
    Vlv_idx2            = find(diff(Vlv_vals_tshort) <= 0.01 & diff(Vlv_vals_tshort) >= -0.01);
    int_idx2            = intersect(Vlv_idx,Vlv_idx2);
    ESP_idx             = int_idx2(find(diff(int_idx2) > 20)+1);
    ESP_vals            = Plv_vals_tshort(ESP_idx);

    % identify appropriate indices for Left Ventricualr End DIastolic Volume
    % (EDV) and End Diastolic Pressure (ESP)
    phi_idx             = find(phi_vals_tshort == 0);
    EDV_idx             = phi_idx(find(diff(phi_idx) > 2));
    EDV_vals            = Vlv_vals_tshort(EDV_idx);

    Vlv_idx3            = find(Vlv_vals_tshort > (meanVlvPks-10));
    int_idx3            = intersect(phi_idx, Vlv_idx3);
    EDP_idx             = int_idx3(find(diff(int_idx3)>2));
    EDP_vals            = Plv_vals_tshort(EDP_idx);

    modelEval.HR        = mean(getsampleusingtime(HR_vals,165,180),'Weighting','time');
    modelEval.CO        = mean(CO_SV.LCO(:,1));
    modelEval.EjFr      = mean(CO_SV.LSV)/mean(Vlv_pks);
    modelEval.SP        = mean(Psys_pks);
    modelEval.DP        = mean(-1*Pdia_pks);
%     modelEval.PP        = mean(Psys_pks - (-1*Pdia_pks));
    modelEval.PP        = mean(Psys_pks) - mean(-1*Pdia_pks);
    modelEval.MAP       = 1/3*mean(Psys_pks) + 2/3*mean(-1*Pdia_pks);
    modelEval.TR        = mean(getsampleusingtime(Rsp_vals,165,180),'Weighting','time') +  mean(getsampleusingtime(Rep_vals,165,180),'Weighting','time') + mean(getsampleusingtime(Rmp_vals,165,180),'Weighting','time');
    modelEval.Emaxlv    = mean(getsampleusingtime(Emaxlv_vals,165,180),'Weighting','time');
    modelEval.ESV       = mean(ESV_vals);
    modelEval.ESP       = mean(ESP_vals);
    modelEval.EDV       = mean(EDV_vals);
    modelEval.EDP       = mean(EDP_vals);
%     modelEval.LSV       = mean(CO_SV.LSV(end-50:end));
    
    %disp('Systolic P, Diastolic P, HR, CO, EjFr')
    simRes(1) = modelEval.SP;
    simRes(2) = modelEval.DP;
    simRes(3) = modelEval.HR; 
    simRes(4) = modelEval.CO;
    simRes(5) = modelEval.EjFr;

    
    %%
    phi_idx             = find(phi_vals == 0);
    EDV_idx             = phi_idx(find(diff(phi_idx) > 2));
    EDV_vals            = Vlv_vals(EDV_idx);
    

    
end