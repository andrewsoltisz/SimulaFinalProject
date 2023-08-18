function [states, varargout] = base_model_init_states_guinea_pig()
  % % Default state values for the guinea pig base model from:
  % % A Tveito, KH Jaeger, MM Maleckar, WR Giles, and S Wall (2020) 
  % % "Computational translation of drug effects from animal experiments to
  % % human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  % --- Default initial state values --- 
  states = zeros(24, 1);

  % --- I_Na ---
  states(1) = 0.00288729380294; % m;
  states(2) = 0.67940212640600; % j;

  % --- I_NaL ---
  states(3) = 0.00034775302411; % mL;
  states(4) = 0.32994883186622; % hL;

  % --- I_Kr ---
  states(5) = 0.00000448646919; % Xr1;
  states(6) = 0.47411587461839; % Xr2;

  % --- I_Ks ---
  states(7) = 0.00352557620095; % x_Ks;

  % --- I_to ---
  states(8) = 0.90835597324176; % q;
  states(9) = 0.00366111004431; % r;

  % --- I_CaL ---
  states(10) = 0.00000232972962; % d;
  states(11) = 0.99587979089285; % f;
  states(12) = 0.00613879236903; % f_Ca_B;

  % --- I_f ---
  states(13) = 0.25148910274039; % xf;

  % --- RyRs ---
  states(14) = 1.0; % r_RyR;

  % --- Ca Concentrations ---
  states(15) = 0.09044528747698; % cn; % network sarcoplasmic reticulum calcium concentration
  states(16) = 0.00003358973872; % cc; % cytosolic calcium concentration
  states(17) = 0.00004095358928; % cd;
  states(18) = 0.00005012089703; % csl;
  states(19) = 0.09021287630541; % cs; % juctional sarcoplasmic reticulum calcium concentration

  % --- Ca Buffer Concentrations ---
  states(20) = 0.00303806018750; % bc;
  states(21) = 0.00489637132038; % bd;
  states(22) = 3.29059736241968; % bs;
  states(23) = 0.02913160218324; % bsl;

  % --- Membrane potential ---
  states(24) = -82.79934405830207; % V_m;

  if nargout == 2

    % --- State names --- 
    state_names = cell(24, 1);

    % --- I_Na ---
    state_names{1} = 'm';
    state_names{2} = 'j';

    % --- I_NaL ---
    state_names{3} = 'mL';
    state_names{4} = 'hL';

    % --- I_Kr ---
    state_names{5} = 'Xr1';
    state_names{6} = 'Xr2';

    % --- I_Ks ---
    state_names{7} = 'x_Ks';

    % --- I_to ---
    state_names{8} = 'q';
    state_names{9} = 'r';

    % --- I_CaL ---
    state_names{10} = 'd';
    state_names{11} = 'f';
    state_names{12} = 'f_Ca_B';

    % --- I_f ---
    state_names{13} = 'xf';

    % --- RyRs ---
    state_names{14} = 'r_RyR';

    % --- Ca Concentrations ---
    state_names{15} = 'ca_nsr'; % network sarcoplasmic reticulum calcium concentration
    state_names{16} = 'ca_i'; % cytosolic calcium concentration
    state_names{17} = 'cd';
    state_names{18} = 'csl';
    state_names{19} = 'ca_jsr'; % junctional sarcoplasmic reticulum calcium concentration

    % --- Ca Buffer Concentrations ---
    state_names{20} = 'bc';
    state_names{21} = 'bd';
    state_names{22} = 'bs';
    state_names{23} = 'bsl';

    % --- Membrane potential ---
    state_names{24} = 'V';
    varargout(1) = {state_names};
  end
end