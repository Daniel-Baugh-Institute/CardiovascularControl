%% James Park 
%% Baroreflex (BR) model based on M. Ursino 1998
%% Shell code to set up parameters for autonomic regulatory controller for circulatory system 
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

global initVals; 

% if exist('BRmodel_ctrls_LBNP.mat') == 2;
%     load BRmodel_ctrls_LBNP;
% else
    % Carotid sinus afferent pathway (i.e. arterial baroreceptors)
    Pn      = 92;         % mmHg
    ka      = 11.758;     % mmHg
    fmin    = 2.52;       % spikes/s
    fmax    = 47.78;      % spikes/s
    tauz    = 6.37;       % s 
    taup    = 2.076;      % s
    
    % Cardiopulmonary afferent pathway
    Ptn     = 10.8;          % mmHg
    Gl      = 3.5;           % spikes/s/s*mmHg
    fmaxl   = 20;            % spikes/s
    fac0    = 10;            % spikes/s (middle value of fac curve)
    kl      = fmaxl/(4*Gl);  % s-mmHg
    taul    = 2;            % s (estimated value) initially guessed 10s
   
    % Lung stretch receptor afferent pathway 
    Vln         = 2.3;        % nominal lung volume
    Vl0         = 1.9;        % initial lung volume
    alpha       = 0.1;        % proportionality relating thoracic P(Pthor) to lung volume
    Gal         = 23.29;      % spikes/s/L
    tau_lung    = 2;          % time constant for lung stretch receptors
    Galh        = 1;
    Galp        = 0.33;       % gain (i.e. weight) of lung stretch receptor contrib. sympathetic effectors
        
    % Sympathetic efferent pathways
    fes_inf   = 2.10;      % spikes/s
    fes0      = 16.11;     % spikes/s
    fesmin    = 2.66;      % spikes/s
    fes_cc    = 0;         % spikes/s
    kes       = 0.0675;    % s
     
    theta_sp  = -4.6;      % spikes/s - offset term for symp. efferents (basal symp. firing rate?)
    
    % Vagal efferent pathway 
    fev0      = 3.2;        % spikes/s
    fev_inf   = 6.3;        % spikes/s
    fev_cc    = 0;          % spikes/s
    fcs0      = 25;         % spikes/s
    kev       = 7.06;       % spikes/s
    
    theta_v  = -0.68;      % or -1.4 spikes/s - offset term for vagal efferents (basal vagal firing rate?)    
    
    
    % Baroreceptor Gains affecting vagal output
    Gab_vag = 1;          % arterial BR gain (weight) to vagal output
    Gac_vag = -1;         % cardiopulmonary BR gain to vagal output
    Gal_vag = -0.103;     % lung stretch receptor gain (i.e. weight) to vagal output
    
    % Afferent Gains to effectors 
    Gabh = 1;             % arterial BR gain on firing frequency to heart
    Gabp = 1;             % arterial BR gain on firing frequency to peripheral resistances
    Gabv = 1;             % arterial BR gain on firing frequency to systemic venous circulation
    Gach = 2;             % cardiopulmonary BR gain on firing frequency to heart
    Gacp = 2.5;           % cardiopulmonary BR gain on firing frequency to peripheral resistances
    Gacv = 0;             % cardiopulmonary BR gain on firing frequency to systemic venous circ.
    
    
    % Effector parameters

    %      [GEmaxlv GEmaxrv  GRsp   GRep   GVusv    GVuev    GTs     GTv    GRmp   GVumv]
    Gx   = [0.475   0.282    0.695  0.653   -265.4  -107.5   -0.13   0.09   2.81   -25];

    %      [tEmmaxlv  tEmaxrv  tRsp  tRep  tVusv  tVuev  tTs  tTv   tRmp  tVumv]
    taux = [8         8        6     6     20     20     2    1.5   6     20];
    
    %      [DEmaxlv  DEmaxrv  DRsp  DRep  DVusv  DVuev  DTs  DTv   DRmp  DVumv]
    Dx   = [2        2        2     2     5      5      2    0.2   2     5]; 
    
    %          [Emaxlv0  Emaxrv0  Rsp0   Rep0   Vusv0   Vuev0  T0    Rmp0  Vumv0]
    initVals = [2.392    1.412    2.49   0.96   1435.4  1247   0.58  4.13  290];

    % OPEN-LOOP simulation parameters
    Pcs_amp    = 0;
    Pcs        = 92;
    Pcs_freq   = 1*2*pi;
    Psa_const  = 92;
    Pra_const  = 0;
    HR_input   = 60; %bpm
    fes_const  = 0;
    fev_const  = 0;
    
    % hemorrhage parameters (blood loss)
    rate    = 0;   % ml/s  (blood volume loss)
    tstart  = 0;   % time to start blood loss 
    ulimit  = 0;   % maximum blood volume loss
 
    save BRmodel_ctrls_LBNP_2001 Pn ka fmin fmax tauz taup ...
        Ptn Gl fmaxl fac0 kl taul...
        Vln Vl0 alpha Gal tau_lung Galp...
        fes_inf fes0 fesmin fes_cc kes theta_sp...
        fev0 fev_inf fev_cc fcs0 kev Gab_vag Gac_vag Gal_vag theta_v... 
        Gabh Gabp Gabv Gach Gacp Gacv Gx Dx taux
% end

load BRmodel_ctrls_LBNP_2001;
    
    
    







    
    