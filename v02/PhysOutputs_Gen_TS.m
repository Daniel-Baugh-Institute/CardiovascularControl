function simRes = PhysOutputs_Gen_TS(simOut)

simRes = zeros(1,6);

if isnan(simOut.time)
    simRes = NaN;
    %modelEval_db = modelEval;
else
    % start & end times for calculating steady state values
    tshort = [165 180];     % chose 15 sec range since this is what is used clinically to determine BPM
    tlong  = [120 180];     % for CO and SV determination
    tplot  = [195 200];     % for plotting purposes

    simTime         = simOut.time;
    tshortIdx       = find(simTime>= tshort(1) & simTime <= tshort(2));
    tlongIdx        = find(simTime>= tlong(1)  & simTime <= tlong(2));
    tplotIdx        = find(simTime>= tplot(1)  & simTime <= tplot(2));


    Emaxlv_vals     = simOut.Emaxlv;
    Tperiod_vals    = simOut.Tperiod;%0.8759*ones(1,length(Emaxlv_vals));% * ones(1,length(Emaxlv_vals)); % edited for constant Tperiod
    HR_vals         = simOut.HR;%1./(Tperiod_vals/60);%* ones(1,length(Emaxlv_vals)); % edited for constant Tperiod
    %Emaxrv_vals     = simRes.signals(4).values;

    Rsp_vals        = simOut.Rsp;
    Rep_vals        = simOut.Rep;
    Rmp_vals        = simOut.Rmp;

%     Vusv_vals       = simRes.signals(8).values;
%     Vuev_vals       = simRes.signals(9).values;
%     Vumv_vals       = simRes.signals(10).values;
%     Vlung_vals      = simRes.signals(11).values;
%     Vla_vals        = simRes.signals(12).values;
    Vlv_vals        = simOut.Vlv;
%     Vra_vals        = simRes.signals(14).values;
%     Vrv_vals        = simRes.signals(15).values;

    Psa_vals        = simOut.Psa;
%     Pmaxlv_vals     = simRes.signals(17).values;
%     Pmaxrv_vals     = simRes.signals(18).values;
%     CPR_input_vals  = simRes.signals(19).values; % Ppv_Pthor

%     Flow_la_vals    = simRes.signals(20).values;
    Flow_lv_vals    = simOut.Flow_lv;
%     Flow_ra_vals    = simRes.signals(22).values;
    Flow_rv_vals    = simOut.Flow_rv;

    Pla_vals        = simOut.Pla;
    Plv_vals        = simOut.Plv;
%     Pra_vals        = simRes.signals(26).values;
%     Prv_vals        = simRes.signals(27).values;

%     BRcellgrp_ff_vals   = simRes.signals(28).values;
%     LScellgrp_ff_vals   = simRes.signals(29).values;
%     CPcellgrp_ff_vals   = simRes.signals(30).values;
% 
%     NAinput_ff_vals     = simRes.signals(31).values;
%     NActr_input_ff_vals = simRes.signals(32).values;
%     DMVinput_ff_vals    = simRes.signals(33).values;

%     fesh_vals       = simRes.signals(34).values;
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

    modelEval.HR        = mean(getsampleusingtime(HR_vals,tshort(1),tshort(2)),'Weighting','time');
    modelEval.CO        = mean(CO_SV.LCO(:,1));
    modelEval.EjFr      = mean(CO_SV.LSV)/mean(Vlv_pks);
    modelEval.SP        = mean(Psys_pks);
    modelEval.DP        = mean(-1*Pdia_pks);
%     modelEval.PP        = mean(Psys_pks - (-1*Pdia_pks));
    modelEval.PP        = mean(Psys_pks) - mean(-1*Pdia_pks);
    modelEval.MAP       = 1/3*mean(Psys_pks) + 2/3*mean(-1*Pdia_pks);
    modelEval.TR        = mean(getsampleusingtime(Rsp_vals,tshort(1),tshort(2)),'Weighting','time') +  mean(getsampleusingtime(Rep_vals,tshort(1),tshort(2)),'Weighting','time') + mean(getsampleusingtime(Rmp_vals,tshort(1),tshort(2)),'Weighting','time');
    modelEval.Emaxlv    = mean(getsampleusingtime(Emaxlv_vals,tshort(1),tshort(2)),'Weighting','time');
    modelEval.ESV       = mean(ESV_vals);
    modelEval.ESP       = mean(ESP_vals);
    modelEval.EDV       = mean(EDV_vals);
    modelEval.EDP       = mean(EDP_vals);
    modelEval.LSV       = mean(CO_SV.LSV(end-50:end));
    
    %disp('Systolic P, Diastolic P, HR, CO, EjFr')
    simRes(1) = modelEval.SP;
    simRes(2) = modelEval.DP;
    simRes(3) = modelEval.HR; %6999999;%69; edited for constant Tperiod
    simRes(4) = modelEval.CO;
    simRes(5) = modelEval.EjFr;
    simRes(6) = modelEval.LSV;
%     disp('Stroke volume for final 15 heartbeats')

% plot(simOut.time(tlongIdx),Psa_vals(tlongIdx))
    
    %%
    phi_idx             = find(phi_vals == 0);
    EDV_idx             = phi_idx(find(diff(phi_idx) > 2));
    EDV_vals            = Vlv_vals(EDV_idx);
    

%     plot(simTime(tshortIdx),Vlv_vals(tshortIdx))
%     hold on
%     plot(simTime(EDV_idx),Vlv_vals(EDV_idx),'o')
%     xlim([simTime(tshortIdx(1)) simTime(tshortIdx(end))])
    
end