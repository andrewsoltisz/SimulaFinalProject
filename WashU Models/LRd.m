function dy=LRd(t,y,p,p_names)

% stimulation settings
st=p(strcmp(p_names,'st'));
fnsh=p(strcmp(p_names,'fnsh'));
Is=p(strcmp(p_names,'Is'));
In=Is*(t>st & t<st+fnsh); 

% extracellular ionic concentrations
ca_o=p(strcmp(p_names,'ca_o'));
na_o=p(strcmp(p_names,'na_o'));
k_o=p(strcmp(p_names,'k_o'));

% physical constants
AF=p(strcmp(p_names,'AF'));
F=p(strcmp(p_names,'F'));
frt=p(strcmp(p_names,'frt'));

% cell geometry
vmyo=p(strcmp(p_names,'vmyo'));
vnsr=p(strcmp(p_names,'vnsr'));
vjsr=p(strcmp(p_names,'vjsr'));

% translation and sensitivity input constants
GNa=p(strcmp(p_names,'g_Na'));
GKsmax=p(strcmp(p_names,'g_Ks'));
GKrmax=p(strcmp(p_names,'g_Kr'));
GK1max=p(strcmp(p_names,'g_K1'));
ibarnak=p(strcmp(p_names,'g_NaK'));
GCaB=p(strcmp(p_names,'g_bCa'));
iupbar=p(strcmp(p_names,'J_SERCA_bar'));
GRel=p(strcmp(p_names,'K_RyR'));
pca=p(strcmp(p_names,'g_CaL'));
vNCX = p(strcmp(p_names,'g_NaCa'));

% other constants
GNab=p(strcmp(p_names,'GNab'));
c1=p(strcmp(p_names,'c1'));
GKpmax=p(strcmp(p_names,'GKpmax'));
ibarpca=p(strcmp(p_names,'ibarpca'));
K_Relss=p(strcmp(p_names,'K_Relss'));
qn=p(strcmp(p_names,'qn'));
tau=p(strcmp(p_names,'tau'));
gacai=p(strcmp(p_names,'gacai'));
gacao=p(strcmp(p_names,'gacao'));
pna=p(strcmp(p_names,'pna'));
ganai=p(strcmp(p_names,'ganai'));
gcat=p(strcmp(p_names,'gcat'));
prnak=p(strcmp(p_names,'prnak'));
kmpca=p(strcmp(p_names,'kmpca'));
nsrbar=p(strcmp(p_names,'nsrbar'));
kmup=p(strcmp(p_names,'kmup'));
gammas=p(strcmp(p_names,'gammas'));
c2=p(strcmp(p_names,'c2'));
kmnai=p(strcmp(p_names,'kmnai'));
kmko=p(strcmp(p_names,'kmko'));
cmdnbar=p(strcmp(p_names,'cmdnbar'));
trpnbar=p(strcmp(p_names,'trpnbar'));
kmtrpn=p(strcmp(p_names,'kmtrpn'));
kmcmdn=p(strcmp(p_names,'kmcmdn'));
csqnbar=p(strcmp(p_names,'csqnbar'));
ganao=p(strcmp(p_names,'ganao'));
pk=p(strcmp(p_names,'pk'));
gaki=p(strcmp(p_names,'gaki'));
gako=p(strcmp(p_names,'gako'));
kmca=p(strcmp(p_names,'kmca'));
kmcsqn=p(strcmp(p_names,'kmcsqn'));
tautr=p(strcmp(p_names,'tautr'));

%Unpack the state vector
V=y(1);
H=y(2);
m=y(3);
J=y(4);
d=y(5);
f=y(6); 
xr=y(7);
ca_T=y(8); 
na_i=y(9);
k_i=y(10);
jsr_T=y(11); 
nsr=y(12);
xs=y(13);
B=y(14);
G=y(15);
xs2=y(16);
Rel=y(17);
Over=y(18);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ca_i]=conc_cai(ca_T,cmdnbar,trpnbar,kmtrpn,kmcmdn);
% [ca_i]=conc_catrpn(ca_T,TRPN,data);

[jsr]=conc_jsr(jsr_T,csqnbar,kmcsqn);

[ina,inab,am,bm,aH,bH,aj,bj] = comp_ina2005(V,m,H,J,na_i,na_o,frt,GNa,GNab);

[ilca,ilcana,ilcak,taud,dss,tauf,fss]=comp_ical2005(V,d,f,ca_i,na_i,k_i,pca,F,frt,gacai,gacao,ca_o,pna,ganai,ganao,na_o,pk,gaki,gako,k_o,kmca);
[inaca]=comp_inaca2000(V,ca_i,na_i,c1,gammas,frt,ca_o,na_o,c2,vNCX);
[inak]=comp_inak2000(V,na_i,ibarnak,na_o,frt,kmnai,kmko,k_o);

[iks,xss,tauxs]=comp_iks2000(V,xs,xs2,ca_i,na_i,k_i,GKsmax,k_o,na_o,prnak,frt);
[icat,bss,gss,taub,taug]=comp_icat2000(V,B,G,ca_i,ca_o,frt,gcat);
[ikr,xrss,tauxr]= comp_ikr95(V,xr,k_i,GKrmax,k_o,frt);
[IK1] = comp_IK194(V,k_i,GK1max,k_o,frt);% time independent IK1
[ikp] = comp_ikp(V,k_i,k_o,GKpmax,frt); % plateau

[ileak,iup,ipca,icab,itr]=calcium_2005(V,nsr,jsr,ca_i,GCaB,ibarpca,ca_o,frt,kmpca,iupbar,nsrbar,kmup,tautr);

caiont = ilca+icab+ipca-2*inaca+icat; %


naiont = ina+inab+3*inaca+ilcana+3*inak;
kiont = ikr+iks+IK1+ikp+ilcak-2*inak-In;%+ito;+insna+insk+ikna

%% Derivatives of state variables first cell

dV= -(naiont+kiont+caiont);

 dH=aH*(1-H)-bH*H;
 dm=am*(1-m)-bm*m; 
 dJ=aj*(1-J)-bj*J;



dD=(dss-d)/taud;
df=(fss-f)/tauf;
dxr=(xrss-xr)/tauxr;
dxs=(xss-xs)/tauxs;
dxs2=(xss-xs2)/tauxs/4;

dnai=-naiont*AF/(vmyo);
dki=-kiont*AF/(vmyo);
dB=(bss-B)/taub;
dG=(gss-G)/taug;




Rel_ss=ilca.*GRel/(1+(K_Relss./jsr).^qn);% alternans with  tau=5 at 300 bcl
tau_Rel=tau./(1+0.0123./jsr);


dRel=-(Rel_ss + Rel)./tau_Rel;



dOver=0;
dcai =-caiont*AF/(vmyo*2)+(ileak-iup)*vnsr/vmyo+(Over+Rel)*vjsr/vmyo;

dnsr = iup-itr*vjsr./vnsr-ileak;%;

djsr = itr-(Rel);

% RETURN DERIVATIVES

 dy = [dV;dH;dm;dJ;dD;df;dxr;dcai;dnai;dki;djsr;dnsr;dxs;dB;dG;dxs2;dRel;dOver];%dCtrap;dcass;dfca;dTRP;dCMD;dLB;dRB;dCSQN];%dyd;dz


%% L-type calcium channel
function [ilca,ilcana,ilcak,taud,dss,tauf,fss]=comp_ical2005(v,d,f,cai,nai,ki,pca,F,frt,gacai,gacao,ca_o,pna,ganai,ganao,na_o,pk,gaki,gako,k_o,kmca)
% Calculates Currents through L-Type Ca Channel

dss=1./(1+exp(-(v+10)/6.24));
taud=dss.*(1-exp(-(v+10)/6.24))./(0.035*(v+10));
dss1=1./(1+exp(-(v+60)/0.024));
dss=dss*dss1;
fss=1./(1+exp((v+32)/8))+(0.6)./(1+exp((50-v)/20));

tauf=1./(0.0197*exp(-(0.0337*(v+10))^2)+0.02);
%ibarca= pca*4*(v*F*frt)*((gacai*cass*exp(2*v*frt)-gacao*ca_o)/(exp(2*v*frt)-1));
ibarca= pca*4*v*F*frt*((gacai*cai*exp(2*v*frt)-gacao*ca_o)/(exp(2*v*frt)-1));

ibarna= pna*(v*F*frt).*((ganai*nai*exp(v*frt)-ganao*na_o)./(exp(v*frt)-1));
ibark= pk*(v*F*frt).*((gaki*ki*exp(v*frt)-gako*k_o)./(exp(v*frt)-1));

fca =1./(1+(cai./kmca));

ilca   = d.*f.*fca.*ibarca;
ilcana = d.*f*fca*ibarna;
ilcak = d.*f*fca*ibark;
%%

%% Fast and background Sodium channels


function [In,inab,am,bm,ah,bh,aj,bj]=comp_ina2005(V,m,H,J,Na_i,na_o,frt,GNa,GNab)


ENa =log(na_o./Na_i)/frt;       % Nernst potential of Na, mV
                    

gNa =GNa*m*m*m*H*J;
In = gNa.*(V-ENa);


inab = GNab*(V-ENa);

a=1-1./(1+exp(-(V+40)/0.024));
ah= a.*0.135.*exp((80+V)./(-6.8));
bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));

aj =  a.*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V)).*(V+37.78)./(1+exp(0.311*(V+79.23)));

bj= (1-a).*(0.3*exp(-2.535e-7*V)./(1+exp(-0.1*(V+32))))+(a).*(0.1212*exp(-0.01052*V)./(1+exp(-0.1378*(V+40.14))));


am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
bm = 0.08*exp(-V/11);


%% Transient calcium channel


function [icat,bss,gss,taub,taug]=comp_icat2000(V,b,g,cai,ca_o,frt,gcat)
%Calculates Currents through T-Type Ca Channel

bss = 1/(1+exp(-(V+14.0)/10.8));
taub = 3.7+6.1/(1+exp((V+25.0)/4.5));
gss = 1/(1+exp((V+60.0)/5.6));

a=1-1./(1+exp(-V/0.0024));
taug = a.*(-0.875*V+12.0)+12.0*(1-a);


ECa = log(ca_o/cai)/2/frt;

icat = gcat*b*b*g*(V-ECa);


%% Time-independent and plato potassium current


function [IK1] = comp_IK194(V,K_i,GK1max,k_o,frt)
% IK1    Time-independent potassium current

EK = log(k_o/K_i)/frt;

ak1 = 1.02/(1+exp(0.2385*(V-EK-59.215)));
bk1 = (0.49124*exp(0.08032*(V-EK+5.476))+exp(0.06175*(V-EK-594.31)))/...
    (1+exp(-0.5143*(V-EK+4.753)));

gK1 = GK1max*sqrt(k_o/5.4)*ak1/(ak1+bk1);
IK1 = gK1*(V-EK);

function [ikp] = comp_ikp(V,K_i,k_o,GKpmax,frt)



EK = log(k_o/K_i)/frt;

ikp = GKpmax*(V-EK)./(1+exp((7.488-V)./5.98)); % plato K 95


%% Slow Activating potassium Current


function [iks,xss,tauxs]=comp_iks2000(v,xs1,xs2,cai,nai,ki,GKsmax,k_o,na_o,prnak,frt)

gks = GKsmax*(1+0.6/(1+(3.8e-5/cai)^1.4));
eks = log((k_o+prnak*na_o)/(ki+prnak*nai))/frt;

xss = 1/(1+exp(-(v-1.5)/16.7));

tauxs = 1/(0.0000719*(v+30)/(1-exp(-0.148*(v+30)))+0.000131*(v+30)/(exp(0.0687*(v+30))-1));


iks = gks*xs1*xs2*(v-eks);
%%
%% Rapidly Activating Potassium Current


function [ikr,xrss,tauxr]= comp_ikr95(v,xr,ki,GKrmax,k_o,frt)

%Calculates Rapidly Activating K Current


gkr = GKrmax*(k_o/5.4).^(1/2);
ekr = log(k_o/ki)/frt;
r = 1/(1+exp((v+9)/22.4));

ikr = gkr*xr*r*(v-ekr);
xrss = 1/(1+exp(-(v+21.5)/7.5));
tauxr = 1/(0.00138*(v+14.2)/(1-exp(-0.123*(v+14.2)))+0.00061*(v+38.9)/(exp(0.145*(v+38.9))-1));
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%  Intracellular Calcium subsystem currents and buffers
function [ileak,iup,ipca,icab,itr]=calcium_2005(v,nsr,jsr,ca_i,GCaB,ibarpca,ca_o,frt,kmpca,iupbar,nsrbar,kmup,tautr)



% NSR Ca Ion Concentration Changes */




ipca = (ibarpca*ca_i)/(kmpca+ca_i);	 % sarcolema pump Ca SERCA
icab =GCaB*(v- log(ca_o/ca_i)/2/frt); % background Ca

ileak = iupbar/nsrbar*nsr;




iup = iupbar.*ca_i/(ca_i+kmup);
itr = (nsr-jsr)./tautr;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculates Na-Ca Exchanger Current

function [Inaca]=comp_inaca2000(v,cai,nai,c1,gammas,frt,ca_o,na_o,c2,vNCX)
%Calculates Na-Ca Exchanger Current


%    inaca;               % NaCa exchanger current (uA/uF)
Inaca = vNCX.*c1*exp(( gammas-1)*v* frt).*((exp(v* frt).*nai.^3*ca_o- na_o^3*cai)./...
    (1+ c2*exp(( gammas-1)*v* frt).*(exp(v* frt).*nai.^3*ca_o+ na_o^3*cai)));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sodium-Potassium Pump
function [inak]=comp_inak2000(v,nai,ibarnak,na_o,frt,kmnai,kmko,k_o)
% Sodium-Potassium Pump */



sigma = (exp(na_o/67.3)-1)/7;

fnak = 1/(1+0.1245*exp((-0.1*v*frt)) + 0.0365*sigma*exp(-v*frt));

inak = ibarnak*fnak./(1+(kmnai./nai).^2)./(1+kmko./k_o);

%%




function [cai]=conc_cai(ca_t,cmdnbar,trpnbar,kmtrpn,kmcmdn)

 	bmyo = cmdnbar+trpnbar-ca_t+kmtrpn+kmcmdn;
	cmyo = kmcmdn*kmtrpn -ca_t*(kmtrpn+kmcmdn)+trpnbar*kmcmdn+cmdnbar*kmtrpn;
	dmyo = -kmtrpn*kmcmdn*ca_t;
%      
     cai =( 2*(bmyo.*bmyo-3*cmyo).^(1/2)/3).*cos(acos((9*bmyo.*cmyo-2*bmyo.*bmyo.*bmyo-27*dmyo)./(2*(bmyo.*bmyo-3*cmyo).^1.5))/3)-(bmyo/3);
function [cajsr]=conc_jsr(ca_t,csqnbar,kmcsqn)

b=csqnbar+kmcsqn-ca_t;
c=ca_t*kmcsqn;
cajsr=-b/2+(b.^2+4*c).^(1/2)/2;

