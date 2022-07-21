%% Define initial values for variables
tstall              = 1;
LBNP_input          = 0; 

% NTS cell groups 
BRcellgrpGainNTS    = 1;
LScellgrpGainNTS    = 1;
CPcellgrpGainNTS    = 1;

BRcellgrpGainDMV    = 0;
LScellgrpGainDMV    = 1;
CPcellgrpGainDMV    = 1;

BR_mdl_params_LBNP_2001
BR_ctrl_params_LBNP_2001

load 'BSparams20160704';
load 'BSparams_LSRest20160715';
load 'BSparams_CPRest20160716';
%load 'BSparams20160704overwrite'; % Edited

BSparams_BRest      = BSparams20160704;
BSparams_LSRest     = BSparams_LSRest20160715;
BSparams_CPRest     = BSparams_CPRest20160716;

% Gains in autonomic system (inverse Emaxlv determination)
Gs_Emaxlv   = BSparams_BRest(25);
Gv_Emaxlv   = BSparams_BRest(26); 

% Psa_input   = 92;

% NTS cell group parameter assignment        variable No.
delta_BR        = BSparams_BRest(1);            %  (1)
minval_BR       = BSparams_BRest(2);            %  (2)
kBR             = BSparams_BRest(3);            %  (3)
midptBR_gain    = BSparams_BRest(4);            %  (4)

delta_LS        = BSparams_LSRest(1);           %  (5)
minval_LS       = BSparams_LSRest(2);           %  (6)
kLS             = BSparams_LSRest(3);           %  (7)
midptLS_gain    = BSparams_LSRest(4);           %  (8)

delta_CP        = BSparams_CPRest(1);           %  (9)
minval_CP       = BSparams_CPRest(2);           % (10)
kCP             = BSparams_CPRest(3);           % (11)
midptCP_gain    = BSparams_CPRest(4);           % (12)

% assign Nucleus Ambiguus parameter assignment
delta_NA        = BSparams_BRest(13);           % (13)
minval_NA       = BSparams_BRest(14);           % (14)
kNA             = BSparams_BRest(15);           % (15)
midptNA_gain    = BSparams_BRest(16);           % (16)

% NA-contractility parameter assignment
delta_NActr     = BSparams_BRest(17);           % (17)
minval_NActr    = BSparams_BRest(18);           % (18)
kNActr          = BSparams_BRest(19);           % (19)
midptNActr_gain = BSparams_BRest(20);           % (20)

% Dorsal Motor Nucleus of Vagus parameter assignment
delta_DMV       = BSparams_BRest(21);           % (21)
minval_DMV      = BSparams_BRest(22);           % (22)
kDMV            = BSparams_BRest(23);           % (23)
midptDMV_gain   = BSparams_BRest(24);           % (24)


% Gains Elasticity & Lung stretch receptors 
Gs_Emaxlv       = BSparams_BRest(25);           % (25)
Gv_Emaxlv       = BSparams_BRest(26);           % (26)

Galh            = BSparams_LSRest(21);          % (27)

% Initial Emaxlv0 
Emaxlv0         = 1.283;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added section

% % NA-Principal Neuron parameter assignment
% delta_NA_PN        = 45.25;           % (37)
% minval_NA_PN       = 4.75;           % (38)
% kNA_PN             = 10.7;           % (39)
% midptNA_PN         = 69;           % (40)
% 
% % DMV-Principal Neuron parameter assignment
% delta_DMV_PN        = 63.6;           % (41)
% minval_DMV_PN       = 1.4;           % (42)
% kDMV_PN             = 15;           % (43)
% midptDMV_PN         = 50;           % (44)
% tau_DMV_PN          = 0.25;           % (45)
% 
% % Local circuit neuron parameter assignment
% delta_LCN        = 5.5;           % (46)
% minval_LCN       = 1.5;           % (47)
% kLCN             = 5;           % (48)
% midptLCN         = 6;           % (49)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need this section for some reason otherwise BRclosedLoopModelEval.m
% returns nothing

% This also prevents ICNmdl_batchrun_tuning.m from re-writing parameters,
% so maybe need to declare simulink global variables?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end added section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%% Define "baseline" set of parameters (based on parameter tuning)
%                      Cla     Vula   Rla      P0lv   kElv    Vulv      Emaxlv0(JP)  kRlv      fes_inf
BaselineLHparams    = [19.23,  25,    2.5e-3,  1.5,    0.014,  16.77,   1.283,       3.75e-4,  2.1];
BaselineBSparams    = [BSparams_BRest(1:4), BSparams_LSRest(1:4), BSparams_CPRest(1:4), BSparams_BRest(13:26), BSparams_LSRest(21)];
% Edited
% BaselineICNparams   = [delta_NA_PN, minval_NA_PN, kNA_PN, midptNA_PN, delta_DMV_PN, ...
%     minval_DMV_PN, kDMV_PN, midptDMV_PN, tau_DMV_PN, delta_LCN, minval_LCN, kLCN, midptLCN];
BaselineNTS_grp_input_idx  = [1 2 3];

%BaselineCVparams = [BaselineBSparams BaselineLHparams BaselineNTS_grp_input_idx];
BaselineCVparams = [BaselineBSparams BaselineLHparams];% BaselineICNparams]; % Edited

%
% gains, time constants, and delays from original Ursino model
%             [GEmaxlv GEmaxrv  GRsp   GRep   GVusv    GVuev    Gs_Tper  Gv_Tper  GRmp   GVumv]
Gx          = [0.475   0.282    0.695  0.653   -265.4  -107.5   -0.13    0.09     2.81   -25];

%             [tEmmaxlv  tEmaxrv  tRsp  tRep  tVusv  tVuev  tTs  tTv   tRmp  tVumv]
taux        = [8         8        6     6     20     20     2    1.5   6     20];

%             [DEmaxlv  DEmaxrv  DRsp  DRep  DVusv  DVuev  DTs  DTv   DRmp  DVumv]
Dx          = [2        2        2     2     5      5      2    0.2   2     5]; 

%             [Emaxlv0  Emaxrv0  Rsp0   Rep0   Vusv0   Vuev0  T0    Rmp0  Vumv0]
initVals    = [2.392    1.412    2.49   0.96   1435.4  1247   0.58  4.13  290];
%}

%% Define global parameters to be manipulated in Rsim 
    
    
% LBNP input value
LBNP_input                  = Simulink.Parameter(LBNP_input);       
LBNP_input.StorageClass     = 'SimulinkGlobal';                     

% NTS cell group parameters
delta_BR                    = Simulink.Parameter(delta_BR);             % (1)
delta_BR.StorageClass       = 'SimulinkGlobal';
minval_BR                   = Simulink.Parameter(minval_BR);            % (2)
minval_BR.StorageClass      = 'SimulinkGlobal';
kBR                         = Simulink.Parameter(kBR);                  % (3)
kBR.StorageClass            = 'SimulinkGlobal';
midptBR_gain                = Simulink.Parameter(midptBR_gain);         % (4)
midptBR_gain.StorageClass   = 'SimulinkGlobal';

delta_LS                    = Simulink.Parameter(delta_LS);             % (5)
delta_LS.StorageClass       = 'SimulinkGlobal';
minval_LS                   = Simulink.Parameter(minval_LS);            % (6)
minval_LS.StorageClass      = 'SimulinkGlobal';
kLS                         = Simulink.Parameter(kLS);                  % (7)
kLS.StorageClass            = 'SimulinkGlobal';
midptLS_gain                = Simulink.Parameter(midptLS_gain);         % (8)
midptLS_gain.StorageClass   = 'SimulinkGlobal';

delta_CP                    = Simulink.Parameter(delta_CP);             % (9)       
delta_CP.StorageClass       = 'SimulinkGlobal';
minval_CP                   = Simulink.Parameter(minval_CP);            % (10)
minval_CP.StorageClass      = 'SimulinkGlobal';
kCP                         = Simulink.Parameter(kCP);                  % (11)
kCP.StorageClass            = 'SimulinkGlobal';
midptCP_gain                = Simulink.Parameter(midptCP_gain);         % (12)
midptCP_gain.StorageClass   = 'SimulinkGlobal';

%%%%%%%%%%%%% ADD THIS SECTION TO INITIALIZEPARAMS.M IF NECESSARY
%%%%%%%%%%%%% %%%%%%%%%%%
% Nucleus Ambiguus
delta_NA                    = Simulink.Parameter(delta_NA);             % (13)
delta_NA.StorageClass       = 'SimulinkGlobal';
minval_NA                   = Simulink.Parameter(minval_NA);            % (14)
minval_NA.StorageClass      = 'SimulinkGlobal';
kNA                         = Simulink.Parameter(kNA);                  % (15)
kNA.StorageClass            = 'SimulinkGlobal';
midptNA_gain                = Simulink.Parameter(midptNA_gain);         % (16)
midptNA_gain.StorageClass   = 'SimulinkGlobal';

% NA-contractility
delta_NActr                     = Simulink.Parameter(delta_NActr);      % (17)
delta_NActr.StorageClass        = 'SimulinkGlobal';
minval_NActr                    = Simulink.Parameter(minval_NActr);     % (18)
minval_NActr.StorageClass       = 'SimulinkGlobal';
kNActr                          = Simulink.Parameter(kNActr);           % (19)
kNActr.StorageClass             = 'SimulinkGlobal';
midptNActr_gain                 = Simulink.Parameter(midptNActr_gain);  % (20)
midptNActr_gain.StorageClass    = 'SimulinkGlobal';

% Dorsal Motor Nucleus of Vagus
delta_DMV                   = Simulink.Parameter(delta_DMV);            % (21)
delta_DMV.StorageClass      = 'SimulinkGlobal';
minval_DMV                  = Simulink.Parameter(minval_DMV);           % (22)
minval_DMV.StorageClass     = 'SimulinkGlobal';
kDMV                        = Simulink.Parameter(kDMV);                 % (23)
kDMV.StorageClass           = 'SimulinkGlobal';
midptDMV_gain               = Simulink.Parameter(midptDMV_gain);        % (24)
midptDMV_gain.StorageClass  = 'SimulinkGlobal';

% Heart elasticity gains for sympathetic and parasympathetic arms
Gs_Emaxlv                   = Simulink.Parameter(Gs_Emaxlv);            % (25)
Gs_Emaxlv.StorageClass      = 'SimulinkGlobal';

Gv_Emaxlv                   = Simulink.Parameter(Gv_Emaxlv);            % (26)
Gv_Emaxlv.StorageClass      = 'SimulinkGlobal';

% Lung Stretch receptor Gain
Galh                        = Simulink.Parameter(Galh);                 % (27)
Galh.StorageClass           = 'SimulinkGlobal';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edited

%     % NA-Principal Neuron parameter assignment
%     delta_NA_PN        = Simulink.Parameter(params2run(37));
%     delta_NA_PN.StorageClass = 'SimulinkGlobal';                % (37)
%     
%     minval_NA_PN       = Simulink.Parameter(params2run(38));    % (38)
%     minval_NA_PN.StorageClass = 'SimulinkGlobal';
%     
%     kNA_PN             = Simulink.Parameter(params2run(39));    % (39)
%     kNA_PN.StorageClass = 'SimulinkGlobal';
%     
%     midptNA_PN         = Simulink.Parameter(params2run(40));    % (40)
%     midptNA_PN.StorageClass = 'SimulinkGlobal';
% 
%     % DMV-Principal Neuron parameter assignment
%     delta_DMV_PN        = Simulink.Parameter(params2run(41));    % (41)
%     delta_DMV_PN.StorageClass = 'SimulinkGlobal';
%     
%     minval_DMV_PN       = Simulink.Parameter(params2run(42));    % (42)
%     minval_DMV_PN.StorageClass = 'SimulinkGlobal';
%     
%     kDMV_PN             = Simulink.Parameter(params2run(43));    % (43)
%     kDMV_PN.StorageClass = 'SimulinkGlobal';
%     
%     midptDMV_PN         = Simulink.Parameter(params2run(44));    % (44)
%     midptDMV_PN.StorageClass = 'SimulinkGlobal';
%     
%     tau_DMV_PN          = Simulink.Parameter(params2run(45));    % (45)
%     tau_DMV_PN.StorageClass = 'SimulinkGlobal';
% 
%     % Local circuit neuron parameter assignment
%     delta_LCN        = Simulink.Parameter(params2run(46));    % (46)
%     delta_LCN.StorageClass = 'SimulinkGlobal';
%     
%     minval_LCN       = Simulink.Parameter(params2run(47));    % (47)
%     minval_LCN.StorageClass = 'SimulinkGlobal';
%     
%     kLCN             = Simulink.Parameter(params2run(48));    % (48)
%     kLCN.StorageClass = 'SimulinkGlobal';
%     
%     midptLCN         = Simulink.Parameter(params2run(49));    % (49)
%     midptLCN.StorageClass = 'SimulinkGlobal';
%     
    % end edited
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameters to modify for testing purposes (recapitulating disease phenotypes)
% IMPORTANT NOTE: If assigning a variable to modulate for testing, need to
% adjust CPR_model4test.m file as well (include as part of rsimsetrtpparam
% construct

% Emaxlv0 input value
Cla                         = Simulink.Parameter(Cla);                  % (28)             
Cla.StorageClass            = 'SimulinkGlobal';

Vula                        = Simulink.Parameter(Vula);                 % (29)
Vula.StorageClass           = 'SimulinkGlobal';

Rla                         = Simulink.Parameter(Rla);                  % (30)
Rla.StorageClass            = 'SimulinkGlobal';

P0lv                        = Simulink.Parameter(P0lv);                 % (31)
P0lv.StorageClass           = 'SimulinkGlobal';

kElv                        = Simulink.Parameter(kElv);                 % (32)
kElv.StorageClass           = 'SimulinkGlobal';

Vulv                        = Simulink.Parameter(Vulv);                 % (33)
Vulv.StorageClass           = 'SimulinkGlobal';

Emaxlv0                     = Simulink.Parameter(Emaxlv0);              % (34)
Emaxlv0.StorageClass        = 'SimulinkGlobal';

kRlv                        = Simulink.Parameter(kRlv);                 % (35)
kRlv.StorageClass           = 'SimulinkGlobal';

fes_inf                     = Simulink.Parameter(fes_inf);               % (36)
fes_inf.StorageClass        = 'SimulinkGlobal';
%}
