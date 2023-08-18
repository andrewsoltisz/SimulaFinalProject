function [values] = base_model_rhs(t, states, parameters, parameterNames)
  % % Compute the right hand side of the human, dog, rabbit and guinea pig 
  % % base models from: A Tveito, KH Jaeger, MM Maleckar, WR Giles, and 
  % % S Wall (2020) "Computational translation of drug effects from animal 
  % % experiments to human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  % Assign states
  if length(states)~=24
    error('Expected the states array to be of size 24.');
  end
  m=states(1); j=states(2); mL=states(3); hL=states(4); Xr1=states(5);...
    Xr2=states(6); x_Ks=states(7); q=states(8); r=states(9); d=states(10);...
    f=states(11); f_Ca_B=states(12); xf=states(13); r_RyR=states(14);...
    cn=states(15); cc=states(16); cd=states(17); csl=states(18);...
    cs=states(19); bc=states(20); bd=states(21); bs=states(22);...
    bsl=states(23); V_m=states(24);

  % Assign parameters
  if length(parameters)~=87
    error('Expected the parameters array to be of size 87.');
  end
  g_Na=parameters(1); lambda_Na=parameters(2); g_NaL=parameters(3);...
    KmKo=parameters(4); KmNaip=parameters(5); g_NaK=parameters(8);...
    epi=parameters(9); g_Ks=parameters(10); pNaK=parameters(11);...
    g_Kr=parameters(14); lambda_K=parameters(15); g_to=parameters(16);...
    g_K1=parameters(17); g_bCl=parameters(18); Q10CaL=parameters(19);...
    g_CaL=parameters(20); Kdact=parameters(21); KmCai=parameters(22);...
    KmCao=parameters(23); KmNai=parameters(24); KmNao=parameters(25);...
    Q10NCX=parameters(26); g_NaCa=parameters(27); ksat=parameters(28);...
    nu=parameters(29); KmPCa=parameters(30); Q10SLCaP=parameters(31);...
    g_pCa=parameters(32); g_bCa=parameters(33); E_f=parameters(34);...
    g_f=parameters(35); Na_i=parameters(36); Na_sl=parameters(37);...
    Nao=parameters(38); K_i=parameters(39); Ko=parameters(40);...
    ce=parameters(41); K_RyR=parameters(42); alpha_RyR=parameters(43);...
    beta_RyR=parameters(44); eta_RyR=parameters(45);...
    gamma_RyR=parameters(46); lambda_RyR=parameters(47); Vc=parameters(48);...
    Vd=parameters(49); Vn=parameters(50); Vs=parameters(51);...
    Vsl=parameters(52); J_SERCA_bar=parameters(53); K_c=parameters(54);...
    K_n=parameters(55); B_tot_c=parameters(56); B_tot_d=parameters(57);...
    B_tot_s=parameters(58); B_tot_sl=parameters(59); k_off_c=parameters(60);...
    k_off_d=parameters(61); k_off_s=parameters(62); k_off_sl=parameters(63);...
    k_on_c=parameters(64); k_on_d=parameters(65); k_on_s=parameters(66);...
    k_on_sl=parameters(67); lambda_B=parameters(68);...
    lambda_B_c=parameters(69); alpha_d_c=parameters(70);...
    alpha_n_s=parameters(71); alpha_sl_c=parameters(72);...
    lambda_c_d=parameters(73); lambda_c_i=parameters(74);...
    lambda_diff=parameters(75); Cli=parameters(76); Clo=parameters(77);...
    Cm=parameters(78); Frdy=parameters(79); R=parameters(80);...
    Temp=parameters(81); chi=parameters(82); lambda_c_e=parameters(83);...
    stim_amplitude=parameters(84); stim_duration=parameters(85);...
    stim_period=parameters(86); stim_start=parameters(87);

  % Init return args
  values = zeros(24, 1);

  % Expressions for the Reversal potentials component
  FoRT = Frdy/(R*Temp);
  ena = log(Nao/Na_i)/FoRT;
  ek = log(Ko/K_i)/FoRT;
  eca_sl = log(ce/csl)/(2*FoRT);
  ecl = log(Cli/Clo)/FoRT;
  Qpow = -31 + Temp/10;

  % Expressions for the I_Na component
  mss = (1 + exp(-19/3 - V_m/9))^(-2);
  taum = 0.06*exp(-(-5/51 + V_m/51)^2) + 0.13*exp(-(23/8 + V_m/16)^2);
  aj = ((V_m >= -40)*(0) + ~(V_m >= -40)*((38 + V_m)*(-25000.0*exp(0.2*V_m) -...
    7e-06*exp(-0.04*V_m))/(1 + 19623624323.7*exp(0.3*V_m))));
  bj = ((V_m >= -40)*(0.6*exp(0.09*V_m)/(1 + exp(-40 - V_m))) + ~(V_m >=...
    -40)*(0.02*exp(-0.01*V_m)/(1 + 0.00369786371648*exp(-0.14*V_m))));
  tauj = 1.0/(aj + bj);
  jss = (1 + exp(72/7 + V_m/7))^(-2);
  values(1) = (-m + mss)/taum;
  values(2) = (-j + jss)/tauj;
  I_Na = g_Na*lambda_Na*m^3*(-ena + V_m)*j;

  % Expressions for the I_NaL component
  mLss = 1.0/(1 + exp(-43/5 - V_m/5));
  tm = 1.0/(8.6*exp(-77/6 - V_m/6) + 6.8*exp(12/35 + V_m/35));
  tmL = tm;
  values(3) = (-mL + mLss)/tmL;
  hLss = 1.0/(1 + 124658.506952*exp(0.133333333333*V_m));
  thL = 200;
  values(4) = (-hL + hLss)/thL;
  GNaL = ((epi == 1)*(g_NaL) + ~(epi == 1)*(0.6*g_NaL));
  I_NaL = lambda_Na*(-ena + V_m)*GNaL*hL*mL;

  % Expressions for the I_NaK component
  sigma = -1/7 + exp(Nao/67)/7;
  fNaK = 1.0/(1 + 0.12*exp(-0.1*FoRT*V_m) + 0.037*exp(-FoRT*V_m)*sigma);
  I_NaK = Ko*g_NaK*fNaK/((1 + KmNaip^4/Na_i^4)*(KmKo + Ko));

  % Expressions for the I_Kr component
  Xr1_inf = 1.0/(1.0 + 0.0146327985189*exp(-0.204081632653*V_m));
  alpha_Xr1 = 450.0/(1.0 + 0.0111089965382*exp(-0.1*V_m));
  beta_Xr1 = 6.0/(1.0 + 13.5813245226*exp(0.0869565217391*V_m));
  tau_Xr1 = 1.0*alpha_Xr1*beta_Xr1;
  values(5) = (-Xr1 + Xr1_inf)/tau_Xr1;
  Xr2_infinity = 1.0/(1.0 + 5.8124373944*exp(0.02*V_m));
  alpha_Xr2 = 3.0/(1.0 + 0.0497870683679*exp(-0.05*V_m));
  beta_Xr2 = 1.12/(1.0 + 0.0497870683679*exp(0.05*V_m));
  tau_Xr2 = 1.0*alpha_Xr2*beta_Xr2;
  values(6) = (-Xr2 + Xr2_infinity)/tau_Xr2;
  I_Kr = 0.430331482912*g_Kr*sqrt(Ko)*(-ek + V_m)*Xr1*Xr2;

  % Expressions for the I_Ks component
  eks = log((Ko + Nao*pNaK)/(K_i + Na_i*pNaK))/FoRT;
  xsss = 1.0/(1 + 0.76228973079*exp(-V_m/14));
  tauxs = 990/(1 + 0.842460441617*exp(-V_m/14));
  values(7) = (-x_Ks + xsss)/tauxs;
  I_Ks = g_Ks*x_Ks^2*(-eks + V_m);

  % Expressions for the I_to component
  q_inf = 1.0/(1.0 + 58.9637634804*exp(0.0769230769231*V_m));
  tau_q = 6 + 39/(0.0168716780457*exp(-0.08*V_m) + 6.46648051673*exp(0.1*V_m));
  values(8) = (-q + q_inf)/tau_q;
  r_inf = 1.0/(1.0 + 3.28489055021*exp(-0.0533333333333*V_m));
  tau_r = 2.75 + 14.4/(0.0207698622486*exp(-0.12*V_m) +...
    15.7194688773*exp(0.09*V_m));
  values(9) = (-r + r_inf)/tau_r;
  I_to = g_to*(-ek + V_m)*q*r;

  % Expressions for the I_K1 component
  aK1 = 1.0/(1 + 7.50455791508e-06*exp(0.2*V_m - 0.2*ek));
  bK1 = (0.745912348821*exp(0.08*V_m - 0.08*ek) +...
    3.32464030033e-16*exp(0.06*V_m - 0.06*ek))/(1 +...
    0.0820849986239*exp(0.5*ek - 0.5*V_m));
  K1ss = aK1/(aK1 + bK1);
  I_K1 = 0.430331482912*g_K1*lambda_K*sqrt(Ko)*(-ek + V_m)*K1ss;

  % Expressions for the I_bCl component
  I_bCl = g_bCl*(-ecl + V_m);

  % Expressions for the I_CaL component
  fss = 1.0/(1 + exp(35/9 + V_m/9)) + 0.6/(1 + exp(5/2 - V_m/20));
  dss = 1.0/(1 + exp(-5/6 - V_m/6));
  taud = (1 - exp(-5/6 - V_m/6))*dss/(0.175 + 0.035*V_m);
  tauf = 1.0/(0.02 + 0.02*exp(-(0.493 + 0.034*V_m)^2));
  values(10) = (-d + dss)/taud;
  values(11) = (-f + fss)/tauf;
  values(12) = -0.012*f_Ca_B + 1.7*(1 - f_Ca_B)*cd;
  ibarca_j = 4*Frdy*g_CaL*(-0.34*ce + 0.34*cd*exp(2*FoRT*V_m))*FoRT*V_m/(-1 +...
    exp(2*FoRT*V_m));
  I_CaL = lambda_c_d*Q10CaL^Qpow*(1 - f_Ca_B)*d*f*ibarca_j;

  % Expressions for the I_NaCa component
  Ka_sl = 1.0/(1 + Kdact^2/csl^2);
  s1_sl = ce*Na_sl^3*exp(nu*FoRT*V_m);
  s2_sl = Nao^3*csl*exp((-1 + nu)*FoRT*V_m);
  s3_sl = KmCao*Na_sl^3 + ce*Na_sl^3 + Nao^3*csl + KmCai*Nao^3*(1 +...
    Na_sl^3/KmNai^3) + KmNao^3*(1 + csl/KmCai)*csl;
  I_NaCa = g_NaCa*lambda_c_e*Q10NCX^Qpow*(-s2_sl + s1_sl)*Ka_sl/((1 +...
    ksat*exp((-1 + nu)*FoRT*V_m))*s3_sl);

  % Expressions for the I_pCa component
  I_pCa = g_pCa*lambda_c_e*Q10SLCaP^Qpow*csl^2/(KmPCa^2 + csl^2);

  % Expressions for the I_bCa component
  I_bCa = g_bCa*lambda_c_e*(-eca_sl + V_m);

  % Expressions for the I_f component
  xf_inf = 1.0/(1 + exp(78/5 + V_m/5));
  tau_xf = 1900/(1 + exp(3/2 + V_m/10));
  values(13) = (-xf + xf_inf)/tau_xf;
  I_f = g_f*(-E_f + V_m)*xf;

  % Expressions for the Ca Fluxes component
  J_CaL = -Cm*chi*I_CaL/(2*Frdy);
  J_pCa = -Cm*chi*I_pCa/(2*Frdy);
  J_bCa = -Cm*chi*I_bCa/(2*Frdy);
  J_NaCa = Cm*chi*I_NaCa/Frdy;
  J_e_sl = J_NaCa + J_bCa + J_pCa;
  J_SERCA = J_SERCA_bar*lambda_c_i*(cc^2/K_c^2 - cn^2/K_n^2)/(1 + cc^2/K_c^2 +...
    cn^2/K_n^2);
  J_n_s = alpha_n_s*lambda_c_i*lambda_diff*(-cs + cn);
  J_sl_c = alpha_sl_c*lambda_c_i*lambda_diff*(-cc + csl);
  J_d_c = alpha_d_c*lambda_c_d*lambda_diff*(-cc + cd);

  % Expressions for the RyRs component
  p = 1.0/(1 + K_RyR^3/cd^3);
  J_RyR_active = alpha_RyR*lambda_RyR*lambda_c_i*(-csl + cs)*p*r_RyR;
  J_leak = alpha_RyR*gamma_RyR*lambda_RyR*lambda_c_i*(-csl + cs);
  values(14) = eta_RyR*(1 - r_RyR)/p -...
    J_RyR_active/(beta_RyR*lambda_RyR*lambda_c_i);
  J_RyR = J_RyR_active + J_leak;

  % Expressions for the Ca Buffers component
  J_c_b = Vc*(-k_off_c*bc + k_on_c*(-bc + B_tot_c*lambda_B*lambda_B_c)*cc);
  J_d_b = Vd*(-k_off_d*bd + k_on_d*(-bd + B_tot_d*lambda_B*lambda_B_c)*cd);
  J_s_b = Vs*(-k_off_s*bs + k_on_s*(-bs + B_tot_s*lambda_B)*cs);
  J_sl_b = Vsl*(-k_off_sl*bsl + k_on_sl*(-bsl +...
    B_tot_sl*lambda_B*lambda_B_c)*csl);

  % Expressions for the Ca Concentrations component
  values(15) = 1.0*(-J_n_s + J_SERCA)/Vn;
  values(16) = 1.0*(-J_SERCA - J_c_b + J_d_c + J_sl_c)/Vc;
  values(17) = 1.0*(-J_d_b - J_d_c + J_CaL)/Vd;
  values(18) = 1.0*(-J_sl_b - J_sl_c + J_RyR + J_e_sl)/Vsl;
  values(19) = 1.0*(-J_RyR - J_s_b + J_n_s)/Vs;

  % Expressions for the Ca Buffer Concentrations component
  values(20) = 1.0*J_c_b/Vc;
  values(21) = 1.0*J_d_b/Vd;
  values(22) = 1.0*J_s_b/Vs;
  values(23) = 1.0*J_sl_b/Vsl;

  % Expressions for the Membrane potential component
  i_Stim = ((V_m < -40)*(1) + ~(V_m < -40)*(0))*((t -...
    stim_period*floor(t/stim_period) <= stim_duration + stim_start & t -...
    stim_period*floor(t/stim_period) >= stim_start)*(-stim_amplitude) + ~(t -...
    stim_period*floor(t/stim_period) <= stim_duration + stim_start & t -...
    stim_period*floor(t/stim_period) >= stim_start)*(0));%
  I_tot = I_CaL + I_K1 + I_Kr + I_Ks + I_Na + I_NaCa + I_NaK + I_NaL + I_bCa...
    + I_bCl + I_f + I_pCa + I_to;
  values(24) = -I_tot - i_Stim;
end
