data.F=96485; %% Faraday constant
data.R=8314;  %% Rey
data.Temp=310;  %% Absolute temperature K
l = 0.01;       % Length of the cell (cm)
a = 0.0011;     % Radius of the cell (cm)

vcell = 1000*pi*a*a*l;     %   3.801e-5 uL   % Cell volume (uL)
ageo = 2*pi*a*a+2*pi*a*l;  %   7.671e-5 cm^2    % Geometric membrane area (cm^2)
Acap = ageo*2;             %   1.534e-4 cm^2    % Capacitive membrane area (cm^2)
data.vmyo = vcell*0.68;    % Myoplasm volume (uL)
vmito = vcell*0.24;  % Mitochondria volume (uL)
data.vsr = vcell*0.06;    % SR volume (uL)
data.vnsr = vcell*0.0552;   % NSR volume (uL)
data.vjsr=  vcell*0.0048;   % JSR volume (uL)
data.vss= vcell*0.02;  %  cell subspace volume (uL)

data.AF=Acap/data.F;

data.frt=data.F/data.Temp/data.R;

data.IKsCa_max=0.6;
data.IKsCa_Kd=38e-6;
 
data.Is=80.0; %%% Stimulus current
data.fnsh=0.5; %% duration of Stimulus current
data.st=0; %% start time of timulus current



data.K_Relss=1;
data.kappa=0.125;
data.tau=4.75;
data.alpha_Rel=data.tau*data.kappa;

data.qn=9;%9
data.tautr=120;%100/12


%%%% SODIUM COMPARTMENT
data.GNa=16;                          % mS/cm^2
data.GNab=0.004; 
data.GNaL=65e-4;

% %% Na-Ca Exchanger Current 2004
data.KmCa=1.25e-4;




data.NCXmax=4.5;
data.ksat=0.27;
data.eta=0.35;
data.KmNai=12.3; data.KmNao=87.5;
data.KmCai=0.0036; data.KmCao=1.3;

% 
% 
% 




%%%% CALCIUM COMPARTMENT


%% L-type channel
data.gacai=1;         % Activity coefficient of Ca
data.gacao=0.341;     % Activity coefficient of Ca

data.kmca=6e-4;     % Half-saturation concentration of Ca channel (mM)
data.pca= 5.4e-4;     % Permiability of membrane to Ca (cm/s)
data.gacai=1;         % Activity coefficient of Ca
data.gacao=0.341;     % Activity coefficient of Ca
data.pna=6.75e-7;     % Permiability of membrane to Na (cm/s)
data.ganai=0.75;      % Activity coefficient of Na
data.ganao=0.75;      % Activity coefficient of Na
data.pk=1.93e-7;       % Permiability of membrane to K (cm/s)
data.gaki=0.75;       % Activity coefficient of K
data.gako=0.75;       % Activity coefficient of K


data.gcat = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SR release




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcium in sarcoplasmic
%%%%%%%%%%%%%%%%%%%%%% reticulum
data.kmup = 0.00092;    % Half-saturation concentration of iup (mM)
data.iupbar = 0.00875;  % Max. current through iup channel (mM/ms)
data.nsrbar = 15;       % Max. [Ca] in NSR (mM)
data.ibarpca = 1.15; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
data.kmpca = 0.5e-3; % Half-saturation concentration of sarcolemmal Ca pump (mM)
data.cmdnbar = 0.050;   % Max. [Ca] buffered in CMDN (mM)
data.trpnbar = 0.070;   % Max. [Ca] buffered in TRPN (mM)
data.kmcmdn = 0.00238;  % Equilibrium constant of buffering for CMDN (mM)
data.kmtrpn = 0.0005;   % Equilibrium constant of buffering for TRPN (mM)
data.trpnf = 40;   % forward  buffered in TRPN (mM)
data.trpnb = 0.02;   % backward  TRPN (mM)
data.cmdnf = 100;   % forward  buffered in TRPN (mM)
data.cmdnb = 0.238;   % backward  TRPN (mM)

data.csqnbar = 10;      % Max. [Ca] buffered in CSQN (mM)
data.kmcsqn = 0.8;      % Equilibrium constant of buffering for CSQN (mM)

data.csqnf = 100;  
data.csqnb = 80; 
data.gcab=0.003016;




data.taudiff=0.2;
data.KmCaMK=0.15;
data.tautr=120;


data.iupmax=0.004375; 
data.Kmup=0.00092;
data.nsrmax=15.0;






data.c1 =0.00025;   % Scaling factor for inaca (uA/uF)
data.c2 = 0.0001;   % Half-saturation concentration of NaCa exhanger (mM)
data.gammas = 0.15;  % Position of energy barrier controlling voltage dependance of inaca






%%%% POTASSIUM COMPARTMENT



data.GKsmax = 0.433;

data.GKrmax =  0.02614;

data.prnak=0.01833;  


data.GK1max=0.75;
data.GKpmax= 0.00552;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sodium-Potassium Pump */



%        inak;    % NaK pump current (uA/uF)
%        fnak;    % Voltage-dependance parameter of inak
%        sigma;   % [Na]o dependance factor of fnak

data.kmnai = 10;    % Half-saturation concentration of NaK pump (mM)
data.kmko = 1.5;    % Half-saturation concentration of NaK pump (mM)
data.ibarnak = 2.25; % Max. current through Na-K pump (uA/uF)

