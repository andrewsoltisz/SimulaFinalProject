function dy=HRd(t,y,p,p_names)

% stimulation settings
start=p(strcmp(p_names,'start')); 
duration=p(strcmp(p_names,'fnsh'));
Is=p(strcmp(p_names,'Is'));

% extracellular ionic concentrations
k_o=p(strcmp(p_names,'Ko'));
na_o=p(strcmp(p_names,'Nao'));
ca_o=p(strcmp(p_names,'Cao'));

% physical constants
F=p(strcmp(p_names,'F'));

% cell geometry
Acap=p(strcmp(p_names,'Acap'));
AF=Acap/F;
vmyo=p(strcmp(p_names,'vmyo'));
vss=p(strcmp(p_names,'vss'));
vnsr=p(strcmp(p_names,'vnsr'));
vjsr=p(strcmp(p_names,'vjsr'));

% translation and sensitivity input constants
%celltype = p(strcmp(p_names,'celltype'));
GK1 = p(strcmp(p_names,'g_K1'));
GKr = p(strcmp(p_names,'g_Kr'));
GKs = p(strcmp(p_names,'g_Ks'));
Gto = p(strcmp(p_names,'g_to'));
GCaL = p(strcmp(p_names,'g_CaL'));
GCab = p(strcmp(p_names,'g_bCa'));
GNa = p(strcmp(p_names,'g_Na'));
GNaL = p(strcmp(p_names,'g_NaL'));
vNKA = p(strcmp(p_names,'g_NaK'));
vNCX = p(strcmp(p_names,'g_NaCa'));
vSERCA = p(strcmp(p_names,'J_SERCA_bar'));
vRyR = p(strcmp(p_names,'K_RyR'));

%Unpack the state vector
V=y(1);
H = y(2);
m = y(3);  
J=y(4); 
d=y(5);
f=y(6); 
xr=y(7); 
ca_i=y(8); 
na_i=y(9);
k_i=y(10);
jsr=y(11);
nsr=y(12);
xs=y(13);
xs2=y(14);
ydv=y(15);
ydv2=y(16);
zdv=y(17);
fca=y(18); 
fca2=y(19); 
f2=y(20);
dp=y(21);
Ctrap=y(22);
ro=y(23);
ri=y(24);
cass=y(25); 
mL=y(26); 
hL=y(27); 
AA=y(28);
CL_i=y(29);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[INa,am,bm,aH,bH,aj,bj] = INatrium2004(V,m,H,J,na_i,na_o,GNa);
[Cactive]=comp_CaMKII(Ctrap,cass);


[ICal,ibarca,dpss,taufca,fcass,taufca2,fca2ss,taud,dss,tauf,fss,tauf2,f2ss]=...
           comp_ical2004(V,d,f,dp,f2,fca,fca2,ca_i,Cactive,cass,na_o,ca_o,GCaL);
[INaCa]=comp_inaca2004(V,ca_i,na_i,na_o,ca_o,vNCX); %%
[INaK]=comp_inak2000(V,na_i,k_o,na_o,vNKA);

[IKs,xss,tauxs]=comp_iks2004(V,xs,xs2,ca_i,na_i,k_i,k_o,na_o,GKs);
[IKr,xrss,tauxr]= comp_ikr2004(V,xr,k_i,k_o,GKr);
[IK1] = IK1_2004(V,k_i,k_o,GK1);
[IKp] = IKp_2004(V,k_i,k_o);

[Ito,ay,by,ay2,by2,ay3,by3]=comp_ito2004(V,k_i,ydv,ydv2,zdv,k_o,Gto);

[INal,amL,bmL,hLss]=comp_inal(V,mL,hL,na_i,na_o,GNaL);
 %% Calcium buffers and time-independent calcium currents
[bmyo,bss,bcsqn,ileak,iup,ipca,ICab,itr,idiff]=calcium_in(V,nsr,jsr,ca_i,cass,Cactive,ca_o,GCab,vSERCA);


[CTNaCl,CTKCl,IClb,Ito2,AAss]=comp_chlor2004(V,AA,CL_i,na_i,k_i,k_o,cass,na_o);

[irelcicr,tauri,riss,ross]=ca_cicr2004(jsr,ICal,ro,ri,ibarca,Cactive,cass,vRyR);

if t>start && t<=start+duration
    In=Is;
else
    In=0;
end

caiont =ICal+ICab+ipca-2*INaCa;%
naiont = INa+3*INaCa+3*INaK+INal;
kiont =IKr+IKs+IK1+IKp-2*INaK+Ito+0.5*In;%+
clont=IClb+Ito2+0.5*In;

%%%% Derivatives of the state variables
dV = -(naiont+kiont+caiont+clont);
%% INa Gates
dH=aH*(1-H)-bH*H;
dm=am*(1-m)-bm*m;
dJ=aj*(1-J)-bj*J;
%% ICaL gates
dD=(dss-d)/taud;
df=(fss-f)/tauf;
df2=(f2ss-f2)/tauf2;
dfca=(fcass-fca)/taufca;
df2ca=(fca2ss-fca2)/taufca2;
dDp=(dpss-dp)/10;%taudp=10
%% IK gates
dxr=(xrss-xr)/tauxr;
dxs=(xss-xs)/tauxs;
dxs2=(xss-xs2)/tauxs/2;
%% Ito gates
dyd=ay*(1-ydv)-by*ydv;
dyd2=ay2*(1-ydv2)-by2*ydv2;
dz=ay3*(1-zdv)-by3*zdv;
%% INaL gates
dmL=amL*(1-mL)-bmL*mL;
dhL=(hLss-hL)/600;
dAA=AAss-AA;%% tauA=1;

dnai=-naiont*AF/(vmyo) + CTNaCl;
dki=-kiont*AF/(vmyo) + CTKCl;
dCLi= clont*AF/(vmyo)+CTNaCl + CTKCl ;  %% valence of chlor is -1

dro=(ross-ro)/3;
dri=(riss-ri)/tauri;

dCtrap=0.05*Cactive*(Cactive-Ctrap)- 6.8e-4*Ctrap;
dcai = bmyo*(-(ICab+ipca-2*INaCa)*AF/(vmyo*2)+(ileak-iup)*vnsr/vmyo+idiff*vss/vmyo);
dcass= bss*(-ICal*AF/(vss*2)+irelcicr*vjsr/vss-idiff);
dnsr = iup-itr*vjsr/vnsr-ileak;
djsr = bcsqn*( itr-irelcicr);

% RETURN DERIVATIVES                   
dy = [dV;dH;dm;dJ;dD;df;dxr;dcai;dnai;dki;djsr;dnsr;dxs;...
    dxs2;dyd;dyd2;dz;dfca;df2ca;df2;dDp;dCtrap;dro;dri;dcass;dmL;dhL;dAA;dCLi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEYOND%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [In,am,bm,ah,bh,aj,bj]=INatrium2004(V,m,H,J,Na_i,na_o,GNa)


frt=0.03743588350780; %FaradayR/T
ENa =log(na_o/Na_i)/frt;       % Nernst potential of Na, mV
%GNa= 8.25;                          % mS/cm^2

gNa = GNa*m*m*m*H*J;
In = gNa*(V-ENa);

if V >= -40
    ah = 0.0;
    aj = 0.0;
    bh = 1/(0.13*(1+exp((V+10.66)/(-11.1))));
    bj = 0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)));
else
    ah = 0.135*exp((80+V)/(-6.8));
    aj = (-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))...
        *(V+37.78)/(1+exp(0.311*(V+79.23)));
    bh = 3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V);
    bj = 0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)));
end

am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
bm = 0.08*exp(-V/11);

function [ical,ibarca,dpss,taufca,fcass,taufca2,fca2ss,taud,dss,tauf,fss,tauf2,f2ss]=comp_ical2004(v,d,f,dp,f2,fca,fca2,cai,Cactive,cass,na_o,ca_o,GCaL)
 frt=0.03743588350780;
F=96485;  
%

%GCaL=2.43e-4;     % Permiability of membrane to Ca (cm/s)
gacai=1;         % Activity coefficient of Ca
gacao=0.341;     % Activity coefficient of Ca
ibarca= GCaL*4*(v-15)*F*frt*((gacai*cass*exp(2*(v-15)*frt)-gacao*ca_o)/(exp(2*(v-15)*frt)-1));
ical= (d.^dp)*f*f2*fca*fca2*ibarca;
dss=1/(1+exp(-(v-4)/6.74));
taud=0.59+0.8*exp(0.052*(v+13))/(1+exp(0.132*(v+13)));

fss=0.7/(1+exp((v+17.12)/7)) + 0.3;
f2ss=0.77/(1+exp((v+17.12)/7)) + 0.23;

tauf= 1/(0.2411*exp(-(0.045*(v-9.6914))^2)+0.0529);
tauf2=1/(0.0423*exp(-(0.059*(v-18.5726))^2)+0.0054);

dpss=9-8/(1+exp(-(v+65)/3.4));


fcass = 0.3/(1-ical/0.05)+0.55/(1+cass/0.003)+0.15;%
fca2ss =1/(1-ical/0.01);

taufca=10*Cactive/(0.15+Cactive)+1/(1+cass/0.003)+0.5;
taufca2=300/(1+exp((-ical-0.175)/0.04))+125;

function IK1 = IK1_2004(V,k_i,k_o,GK1)
% IK1     Time-independent potassium current

frt=0.03743588350780;

GK1 = GK1*sqrt(k_o/5.4);
EK1 = log(k_o/k_i)/frt;

ak1 = 1.02/(1+exp(0.2385*(V-EK1-59.215)));
bk1 = (0.49124*exp(0.08032*(V-EK1+5.476))+exp(0.06175*(V-EK1-594.31)))/...
    (1+exp(-0.5143*(V-EK1+4.753)));

gK1 = GK1*ak1/(ak1+bk1);
IK1 = gK1*(V-EK1);

function IKp=IKp_2004(V,k_i,k_o)

GKp=0.00276;
frt=0.03743588350780;
EK1 = log(k_o/k_i)/frt;
IKp = GKp*(V-EK1)/(1+exp((7.488-V)/5.98)); % K 2004

function [iks,xss,tauxs]=comp_iks2004(v,xs1,xs2,cai,nai,ki,k_o,na_o,GKs)
%Calculates Slowly Activating K Current

prnak=0.01833;  frt=0.03743588350780;

%GKs = 0.0248975;
gks = GKs*(1+0.6/(1+(3.8e-5/cai)^(1.4)));
eks = log((k_o+prnak*na_o)/(ki+prnak*nai))/frt;

xss = 1/(1+exp(-(v-10.5)/24.7));

tauxs = 1/(7.61e-5*(v+44.6)/(1-exp(-9.97*(v+44.6)))+3.6e-4*(v-0.55)/(exp(0.128*(v-0.55))-1));


iks = gks*xs1*xs2*(v-eks);


function [ikr,xrss,tauxr]= comp_ikr2004(v,xr,ki,k_o,GKr)
%Calculates Rapidly Activating K Current
frt=0.03743588350780;

%GKr = 0.0138542;
gkr = GKr*sqrt(k_o/5.4);
ekr = log(k_o/ki)/frt;
r = 1/(1+exp((v+10)/15.4));

ikr = gkr*xr*r*(v-ekr);
xrss = 1/(1+exp(-(v+10.085)/4.25));
tauxr = 1/(6e-4*(v-1.7384)/(1-exp(-0.136*(v-1.7384)))...
    +3e-4*(v+38.3608)/(exp(0.1522*(v+38.3608))-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bmyo,bss,bcsqn,ileak,iup,ipca,icab,itr,idiff]=...
            calcium_in(v,nsr,jsr,ca_i,cass,Cactive,ca_o,GCaB,vSERCA)


F=96485;        % Faraday's Constant (C/mol)


frt=0.03743588350780;

% Max. [Ca] in NSR (mM)
vPMCA =0.0575; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
kmpca = 0.5e-3; % Half-saturation concentration of sarcolemmal Ca pump (mM)
ipca = vPMCA*ca_i/(kmpca+ca_i);	 % sarcolema pump Ca SERCA

kmt=0.5e-3;    kmc=2.38e-3;    tbar=70e-3;    cbar=50e-3;
kmcsqn=0.8;csqnbar=10;
bcsqn = 1/(1+kmcsqn*csqnbar/(jsr+kmcsqn)^2);

bmyo=1/(1+ cbar*kmc/(ca_i+kmc)^2+kmt*tbar/(ca_i+kmt)^2);

% background ICa
%GCaB=1.995084e-7;
icab=GCaB*4*v*F*frt*(ca_i*exp(2*v*frt)-0.341*ca_o)/(exp(2*v*frt)-1);



idiff=(cass-ca_i)/0.2;
KmCaMK=0.15;
itr = (nsr-jsr)/120;
dKmPLBmax=0.00017;
dJupmax=0.75;
dKmPLB=dKmPLBmax*Cactive/(KmCaMK+Cactive);
dJup=dJupmax*Cactive/(KmCaMK+Cactive);
%vSERCA=0.004375; 
Kmup=0.00092; nsrmax=15.0;

iup=(dJup+1.0)*vSERCA*ca_i/(ca_i+Kmup-dKmPLB);

ileak = vSERCA*nsr/15;

BSRmax=0.047; KmBSR=0.00087; 
BSLmax=1.124; KmBSL=0.0087;
bss=1/(1+BSRmax*KmBSR/(KmBSR+cass)^2+BSLmax*KmBSL/(KmBSL+cass)^2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [irelcicr,tauri,riss,ross]=ca_cicr2004(jsr,ical,ro,ri,ibarca,Cactive,Cass,vRyR)
cafac=1/(1+exp((ical+0.05)/0.015));

vg=1/(1+exp((ibarca+13)/5));
KmCaMK=0.15;
%vRyR = 3000;
Grel=3000*vRyR;
dro_inf=(jsr^1.9)/(jsr.^1.9+(49.28*Cass/(Cass+0.0028)).^1.9);

dtau_rel_max=10;
dtau_rel=10*Cactive/(KmCaMK+Cactive);


ross=dro_inf/(1/ical/ical+1);

riss=1/(1+exp((Cass- (4e-4)+0.002*cafac)/2.5e-5));
tauri=3+dtau_rel+(350.0-dtau_rel)/(1.0+exp((Cass-0.003+0.003*cafac)/2e-4));

irelcicr=Grel*ro*ri*(jsr-Cass);




function [ INaCa]=comp_inaca2004(v,cai,nai,na_o,ca_o,vNCX)
%Calculates Na-Ca Exchanger Current
KmCa=1.25e-4;
cai=1.5*cai;
frt=0.03743588350780;
allo=1/(1+(KmCa/cai).^2);




%vNCX=4.5;
ksat=0.27;
eta=0.35;
KmNai=12.3; KmNao=87.5;
KmCai=0.0036; KmCao=1.3;


num=vNCX*(nai^3*ca_o*exp(eta*v*frt)-na_o^3*cai*exp((eta-1)*v*frt));
denom1=1+ksat*exp((eta-1)*v*frt);
denom2=KmCao*nai^3+KmNao^3*cai+KmNai^3*ca_o*(1+cai/KmCai);
denom3=KmCai*na_o^3*(1+(nai/KmNai)^3)+nai^3*ca_o+na_o^3*cai;
dE=num/(denom1*(denom2+denom3));

INaCa=allo*dE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inak]=comp_inak2000(v,nai,k_o,na_o,vNKA)
% Sodium-Potassium Pump */

frt=0.03743588350780;
% Faraday's Constant (C/mol)

%        inak;    % NaK pump current (uA/uF)
%        fnak;    % Voltage-dependance parameter of inak
%        sigma;   % [Na]o dependance factor of fnak

kmnai = 10;    % Half-saturation concentration of NaK pump (mM)
kmko = 1.5;    % Half-saturation concentration of NaK pump (mM)
%vNKA = 0.61875; % Max. current through Na-K pump (uA/uF)

sigma = (exp(na_o/67.3)-1)/7;

fnak = 1/(1+0.1245*exp((-0.1*v*frt)) + 0.0365*sigma*exp(-v*frt));

inak = vNKA*fnak*(1/(1+(kmnai/nai)^2))*(k_o/(k_o+kmko));



%

function [ito,a,b,ai,bi,ai2,bi2]=comp_ito2004(v,ki,y1,y2,z,k_o,Gto)

frt=0.03743588350780;
%Gto = 0.19;
ekdv = log(k_o/ki)/frt;
rv = exp(v/300);

a = 25*exp((v-40)/25)/(1+exp((v-40)/25));
b = 25*exp(-(v+90)/25)/(1+exp(-(v+90)/25));

ai = 0.03/(1+exp((v+60)/5));
bi = 0.2*exp((v+25)/5)/(1+exp((v+25)/5));

ai2 = 0.00225/(1+exp((v+60)/5));
bi2 = 0.1*exp((v+25)/5)/(1+exp((v+25)/5));

ito=Gto*y1^3*y2*z*rv*(v-ekdv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cactive]=comp_CaMKII(Ctrap,cass)

CaMK0=0.05;
Km=0.0015;
Cbound=CaMK0*(1-Ctrap)/(1+Km/cass);
Cactive=Cbound+Ctrap;


%%%%%%%%%%%%%%%%%%%%%%%%%%  CHLOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CTNaCl,CTKCl,IClb,Ito2,AAss]=comp_chlor2004(v,AA,Cli,Nai,Ki,Ko,Cass,na_o)

frt=0.03743588350780;F=96485;
Clo=100;
Kmto2=0.1502;
GClB=2.25e-4;
ENa=log(na_o/Nai)/frt;
EK=log(Ko/Ki)/frt;

ECl=-log(Clo/Cli)/frt;
CTKClmax=7.0756e-6;
CTNaClmax=9.8443e-6;
CTKCl=CTKClmax*(EK-ECl)/((EK-ECl)+87.8251);
%CTKCl=CTKClmax./(1+87.8251./(EK-ECl));
CTNaCl=CTNaClmax*(ENa-ECl)^4.0/((ENa-ECl)^4.0+87.8251^4.0);

%CTNaCl=CTNaClmax./(1+(87.8251./(ENa-ECl)).^4.0);

PCl=4e-7;
Ito2_max=PCl*v*F*frt*(Cli-Clo*exp(v*frt))/(1-exp(v*frt));
AAss=1/(1.0+Kmto2/Cass);

Ito2=Ito2_max*AA;
IClb=GClB*(v-ECl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inal,amL,bmL,hLss]=comp_inal(v,mL,hL,Nai,na_o,GNaL)

frt=0.03743588350780;
%GNaL=65e-4;
ENa=log(na_o/Nai)/frt;
inal=GNaL*mL^3*hL*(v-ENa);
amL=0.32*(v+47.13)/(1-exp(-0.1*(v+47.13)));
bmL=0.08*exp(-v/11.0);

hLss=1/(1+exp((v+91)/6.1));
