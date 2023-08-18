function [parameters, varargout] = base_model_init_parameters_rabbit()
% % Default parameter values for the rabbit base model from:
  % % A Tveito, KH Jaeger, MM Maleckar, WR Giles, and S Wall (2020) 
  % % "Computational translation of drug effects from animal experiments to
  % % human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(87, 1);

  % --- I_Na ---
  parameters(1) = 2.26955216457622; % g_Na;
  parameters(2) = 1; % lambda_Na;

  % --- I_NaL ---
  parameters(3) = 0.02685306332202; % g_NaL;

  % --- I_NaK ---
  parameters(4) = 1.5; % KmKo;
  parameters(5) = 11; % KmNaip;
  parameters(6) = 2; % q_SERCA;
  parameters(7) = 1.6; % Q10NaK;
  parameters(8) = 2.77779344601704; % g_NaK;

  % --- I_Ks ---
  parameters(9) = 1; % epi;
  parameters(10) = 0.01563106875281; % g_Ks;
  parameters(11) = 0.018; % pNaK;

  % --- I_Kr ---
  parameters(12) = 0.025; % L0;
  parameters(13) = 2.3; % Q;
  parameters(14) = 0.08738856250515; % g_Kr;
  parameters(15) = 1; % lambda_K;

  % --- I_to ---
  parameters(16) = 0.17169127018627; % g_to;

  % --- I_K1 ---
  parameters(17) = 0.35418180283185; % g_K1;

  % --- I_bCl ---
  parameters(18) = 0.01103888114719; % g_bCl;

  % --- I_CaL ---
  parameters(19) = 1.8; % Q10CaL;
  parameters(20) = 0.25049604967787; % g_CaL;

  % --- I_NaCa ---
  parameters(21) = 0.00015; % Kdact;
  parameters(22) = 0.0036; % KmCai;
  parameters(23) = 1.3; % KmCao;
  parameters(24) = 12.3; % KmNai;
  parameters(25) = 87.5; % KmNao;
  parameters(26) = 1.6; % Q10NCX;
  parameters(27) = 9.78290640505123; % g_NaCa;
  parameters(28) = 0.3; % ksat;
  parameters(29) = 0.3; % nu;

  % --- I_pCa ---
  parameters(30) = 0.0005; % KmPCa;
  parameters(31) = 2.35; % Q10SLCaP;
  parameters(32) = 0.06444216194260; % g_pCa;

  % --- I_bCa ---
  parameters(33) = 0.00050436712647; % g_bCa;

  % --- I_f ---
  parameters(34) = -17; % E_f;
  parameters(35) = 0.0001; % g_f;

  % --- Na Concentrations ---
  parameters(36) = 8; % Na_i;
  parameters(37) = 8; % Na_sl;
  parameters(38) = 140; % Nao;

  % --- K Concentration ---
  parameters(39) = 120.0; % K_i;
  parameters(40) = 5.4; % Ko;

  % --- Ca Concentrations ---
  parameters(41) = 1.8; % ce;

  % --- RyRs ---
  parameters(42) = 0.015; % K_RyR;
  parameters(43) = 0.0075; % alpha_RyR;
  parameters(44) = 0.038; % beta_RyR;
  parameters(45) = 1e-05; % eta_RyR;
  parameters(46) = 0.001; % gamma_RyR;
  parameters(47) = 1; % lambda_RyR;

  % --- Intracellular volumes ---
  parameters(48) = 0.917; % Vc;
  parameters(49) = 0.001; % Vd;
  parameters(50) = 0.05; % Vn;
  parameters(51) = 0.004; % Vs;
  parameters(52) = 0.028; % Vsl;

  % --- SERCA pump ---
  parameters(53) = 0.00024; % J_SERCA_bar;
  parameters(54) = 0.00025; % K_c;
  parameters(55) = 1.7; % K_n;

  % --- Ca Buffers ---
  parameters(56) = 0.07; % B_tot_c;
  parameters(57) = 1.2; % B_tot_d;
  parameters(58) = 27; % B_tot_s;
  parameters(59) = 0.9; % B_tot_sl;
  parameters(60) = 0.03; % k_off_c;
  parameters(61) = 1; % k_off_d;
  parameters(62) = 65.0; % k_off_s;
  parameters(63) = 0.15; % k_off_sl;
  parameters(64) = 40; % k_on_c;
  parameters(65) = 100; % k_on_d;
  parameters(66) = 100.0; % k_on_s;
  parameters(67) = 100; % k_on_sl;
  parameters(68) = 1; % lambda_B;
  parameters(69) = 1; % lambda_B_c;

  % --- Ca Fluxes ---
  parameters(70) = 0.0017; % alpha_d_c;
  parameters(71) = 0.012; % alpha_n_s;
  parameters(72) = 0.15; % alpha_sl_c;
  parameters(73) = 1; % lambda_c_d;
  parameters(74) = 1; % lambda_c_i;
  parameters(75) = 1; % lambda_diff;

  % --- Cl Concentrations ---
  parameters(76) = 15; % Cli;
  parameters(77) = 150; % Clo;

  % --- Membrane potential ---
  parameters(78) = 0.01; % Cm;
  parameters(79) = 96.485; % Frdy;
  parameters(80) = 8.314; % R;
  parameters(81) = 310; % Temp;
  parameters(82) = 0.6; % chi;
  parameters(83) = 1; % lambda_c_e;
  parameters(84) = 40.0; % stim_amplitude;
  parameters(85) = 20.0; % stim_duration;
  parameters(86) = 1000.0; % stim_period;
  parameters(87) = 50.0; % stim_start;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(87, 1);

    % --- I_Na ---
    parameter_names{1} = 'g_Na';
    parameter_names{2} = 'lambda_Na';

    % --- I_NaL ---
    parameter_names{3} = 'g_NaL';

    % --- I_NaK ---
    parameter_names{4} = 'KmKo';
    parameter_names{5} = 'KmNaip';
    parameter_names{6} = 'q_SERCA';
    parameter_names{7} = 'Q10NaK';
    parameter_names{8} = 'g_NaK';

    % --- I_Ks ---
    parameter_names{9} = 'epi';
    parameter_names{10} = 'g_Ks';
    parameter_names{11} = 'pNaK';

    % --- I_Kr ---
    parameter_names{12} = 'L0';
    parameter_names{13} = 'Q';
    parameter_names{14} = 'g_Kr';
    parameter_names{15} = 'lambda_K';

    % --- I_to ---
    parameter_names{16} = 'g_to';

    % --- I_K1 ---
    parameter_names{17} = 'g_K1';

    % --- I_bCl ---
    parameter_names{18} = 'g_bCl';

    % --- I_CaL ---
    parameter_names{19} = 'Q10CaL';
    parameter_names{20} = 'g_CaL';

    % --- I_NaCa ---
    parameter_names{21} = 'Kdact';
    parameter_names{22} = 'KmCai';
    parameter_names{23} = 'KmCao';
    parameter_names{24} = 'KmNai';
    parameter_names{25} = 'KmNao';
    parameter_names{26} = 'Q10NCX';
    parameter_names{27} = 'g_NaCa';
    parameter_names{28} = 'ksat';
    parameter_names{29} = 'nu';

    % --- I_pCa ---
    parameter_names{30} = 'KmPCa';
    parameter_names{31} = 'Q10SLCaP';
    parameter_names{32} = 'g_pCa';

    % --- I_bCa ---
    parameter_names{33} = 'g_bCa';

    % --- I_f ---
    parameter_names{34} = 'E_f';
    parameter_names{35} = 'g_f';

    % --- Na Concentrations ---
    parameter_names{36} = 'Na_i';
    parameter_names{37} = 'Na_sl';
    parameter_names{38} = 'Nao';

    % --- K Concentration ---
    parameter_names{39} = 'K_i';
    parameter_names{40} = 'Ko';

    % --- Ca Concentrations ---
    parameter_names{41} = 'ce';

    % --- RyRs ---
    parameter_names{42} = 'K_RyR';
    parameter_names{43} = 'alpha_RyR';
    parameter_names{44} = 'beta_RyR';
    parameter_names{45} = 'eta_RyR';
    parameter_names{46} = 'gamma_RyR';
    parameter_names{47} = 'lambda_RyR';

    % --- Intracellular volumes ---
    parameter_names{48} = 'Vc';
    parameter_names{49} = 'Vd';
    parameter_names{50} = 'Vn';
    parameter_names{51} = 'Vs';
    parameter_names{52} = 'Vsl';

    % --- SERCA pump ---
    parameter_names{53} = 'J_SERCA_bar';
    parameter_names{54} = 'K_c';
    parameter_names{55} = 'K_n';

    % --- Ca Buffers ---
    parameter_names{56} = 'B_tot_c';
    parameter_names{57} = 'B_tot_d';
    parameter_names{58} = 'B_tot_s';
    parameter_names{59} = 'B_tot_sl';
    parameter_names{60} = 'k_off_c';
    parameter_names{61} = 'k_off_d';
    parameter_names{62} = 'k_off_s';
    parameter_names{63} = 'k_off_sl';
    parameter_names{64} = 'k_on_c';
    parameter_names{65} = 'k_on_d';
    parameter_names{66} = 'k_on_s';
    parameter_names{67} = 'k_on_sl';
    parameter_names{68} = 'lambda_B';
    parameter_names{69} = 'lambda_B_c';

    % --- Ca Fluxes ---
    parameter_names{70} = 'alpha_d_c';
    parameter_names{71} = 'alpha_n_s';
    parameter_names{72} = 'alpha_sl_c';
    parameter_names{73} = 'lambda_c_d';
    parameter_names{74} = 'lambda_c_i';
    parameter_names{75} = 'lambda_diff';

    % --- Cl Concentrations ---
    parameter_names{76} = 'Cli';
    parameter_names{77} = 'Clo';

    % --- Membrane potential ---
    parameter_names{78} = 'Cm';
    parameter_names{79} = 'Frdy';
    parameter_names{80} = 'R';
    parameter_names{81} = 'Temp';
    parameter_names{82} = 'chi';
    parameter_names{83} = 'lambda_c_e';
    parameter_names{84} = 'stim_amplitude';
    parameter_names{85} = 'stim_duration';
    parameter_names{86} = 'stim_period';
    parameter_names{87} = 'stim_start';
    varargout(1) = {parameter_names};
  end
end
