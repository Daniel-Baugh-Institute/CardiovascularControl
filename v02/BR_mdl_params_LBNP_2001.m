%% James Park Feb 2016
%% Baroreflex (BR) model based on M. Ursino 1998
%% Shell code to set up parameters for circulatory system 
%% implementation based on Zhe Hu (HSP model which uses M. Ursino 1998)

%%=====================================
% subscript meanings:
% sa  = systemic arteries
% sp  = splanchic periphery (i.e. abdominal areas)
% ep  = extrasplanchic peripheral circulation (i.e. non abdominal areas)
% ev  = extrasplanchic venous circulation
% ra  = right atrium
% rv  = right ventricle
% la  = left atrium
% lv  = left ventricle
% pa  = pulmonary arteries
% pv  = pulmonary veins
% pp  = pulmonary peripheral circulation
% or  = output from right ventricle
% ol  = output from left ventricle.
% uxx = unstressed volume of compartment xx (e.g. ulv = unstressed volume of left ventricle
% es  = efferent sympathetic pathway
% ev  = efferent vagal pathways
% cs  = carotid sinus (i.e. from baroreceptors)
%%=====================================


%%======================================
% Initialize parameters
%%======================================

global initVals;      % initial vals for baroreflex (BR) model

% if exist('BRmodel_params_LBNP.mat')== 2; 
%     load BRmodel_params_LBNP.mat;
% else
    % Compliance values for circulatory compartments [mL/mmHg]
    Csa    = 0.28;
    Csp    = 2.05;
    Cep    = 1.36;        % based on Ursino 2001 (univentricular flow model) 
%     Cep    = 1.67; 
%     Csv    = 61.11;
%     Cev    = 50.0;
    Csv    = 43.11;       % based on Ursino 2001 (univentricular flow model)    
    Cev    = 28.4;        % based on Ursino 2001 (univentricular flow model)
    Ctv    = 33;          % based on Ursino 2001 (univentricular flow model) 
    Cmp    = 0.31;
    Cmv    = 6.6;     
    Cpa    = 0.76;
    Cpp    = 5.8;
    Cpv    = 25.37;
    
    % Unstressed volumes for circulatory compartments [mL]
    Vusa   = 0;
    Vusp   = 274.4;
    Vuep   = 274.1;      % 336.6;
    Vump   = 62.5; 
    Vusv   = 1121;
    Vuev   = 1120;      % 1375;
    Vumv   = 255;   
    Vupa   = 0;
    Vupp   = 123;
    Vupv   = 120;
    Vutv   = 0; 
    
    % Hydraulic resistance  [mmHg-mL/s]
    Rsa    = 0.06;
    Rsp    = 3.307;
    Rep    = 1.725;       % 1.407;
    Rev    = 0.0197;      % 0.016;
    Rmv    = 0.0848;      
    Rsv    = 0.038;
    Rtv    = 0.0054;      % based on Ursino 2001 (univentricular flow model)
    Rpa    = 0.023;
    Rpp    = 0.0894;
    Rpv    = 0.0056;
    Rd     = 1e4;  
    
    % Inertance  [mmHg-mL/s^2]
    Lsa    = 0.22e-3;
    Lpa    = 0.18e-3;
    
    % Left heart parameters
    Cla     = 19.23;      % mL/mmHg   - compliance
    Vula    = 25;         % mL/mmHg   - unstressed left atrium volume
    Rla     = 2.50e-3;    % mmHg-s/mL - resistance 
    P0lv    = 1.5;        % mmHg-s/mL - basal left ventricle pressure
    kElv    = 0.014;      % 1/mL      - monoexponential constant of pressure during diastole
    Vulv    = 16.77;      % mL        - unstressed left ventricle volume
    Emaxlv  = 2.95;       % mmHg/mL   - elastance left ventricle
    kRlv    = 3.75e-4;     % s/mL     - rate constant
    
    % Right heart parameters
    Cra     = 31.25;      % mL/mmHg   - compliance
    Vura    = 25;         % mL/mmHg   - unstressed left atrium volume
    Rra     = 2.5e-3;     % mmHg-s/mL - resistance
    P0rv    = 1.5;        % mmHg-s/mL - basal right ventricle pressure
    kErv    = 0.011;      % 1/mL      - monoexponential constant of pressure during diastole
    Vurv    = 40.8;       % mL        - unstressed right ventricle volume
    Emaxrv  = 1.75;       % mmHg/mL   - elastance right ventricle 
    kRrv    = 1.4e-3;     % s/mL      - rate constant
    
    % Total volume of blood
    Vt      = 5300;       % mL 
    
    % Heart period/beat
%     Tsys0   = 0.5;        % s
    Tsys0   = 0.4;
    ksys    = 0.075;      
    
    % Respiration parameters
    InspResp = 1.6/1.4;   % ratio of inspiration/respiration 
    InspFrac = 1.6/4.0;   % inspiration to total respiration cycle
    RestFrac = 1/4;       % rest to total respiration lenght
    
    Tresp    = 4.0; 
    Tinsp    = Tresp*InspFrac;
    Texp     = Tinsp/InspResp;  
    
    PthorMax = -4;        % mmHg
    PthorMin = -9;        % mmHg
    
    % Lower body negative pressure
    LBNP    = 0;          % initial lower body negative pressure
    
    save BRmodel_params_LBNP_2001 Csa Csp Cep Csv Cev Ctv Cmp Cmv Cpa Cpp Cpv...
        Vusa Vusp Vuep Vump Vusv Vuev Vumv Vupa Vupp Vupv Vutv...
        Rsa Rsp Rep Rev Rmv Rsv Rtv Rpa Rpp Rpv Rd Lsa Lpa ...
        Cla Vula Rla P0lv kElv Vulv Emaxlv kRlv...
        Cra Vura Rra P0rv kErv Vurv Emaxrv kRrv Vt...
        Tsys0 ksys InspResp InspFrac RestFrac...
        Tresp Tinsp Texp PthorMax PthorMin

% end


load BRmodel_params_LBNP_2001;


