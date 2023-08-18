function [states, varargout] = base_model_init_states_human()
  % % Default state values for the human base model from:
  % % A Tveito, KH Jaeger, MM Maleckar, WR Giles, and S Wall (2020) 
  % % "Computational translation of drug effects from animal experiments to
  % % human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  % --- Default initial state values --- 
  states = zeros(24, 1);

  % --- I_Na ---
  states(1) = 0.00565566388187; % m;
  states(2) = 0.55814710545729; % j;

  % --- I_NaL ---
  states(3) = 0.00066380502723; % mL;
  states(4) = 0.23400507851859; % hL;

  % --- I_Kr ---
  states(5) = 0.00054591339540; % Xr1;
  states(6) = 0.45801990430307; % Xr2;

  % --- I_Ks ---
  states(7) = 0.00443791616325; % x_Ks;

  % --- I_to ---
  states(8) = 0.88542405924077; % q;
  states(9) = 0.00434745567832; % r;

  % --- I_CaL ---
  states(10) = 0.00000399402304; % d;
  states(11) = 0.99390430207136; % f;
  states(12) = 0.01481547054494; % f_Ca_B;

  % --- I_f ---
  states(13) = 0.16493522807906; % xf;

  % --- RyRs ---
  states(14) = 0.99999877286213; % r_RyR;

  % --- Ca Concentrations ---
  states(15) = 0.58928367645190; % cn; % network sarcoplasmic reticulum calcium concentration
  states(16) = 0.00009317813167; % cc; % cytosolic calcium concentration
  states(17) = 0.00010238636602; % cd;
  states(18) = 0.00011412748875; % csl;
  states(19) = 0.58895033791647; % cs; % junctional sarcoplasmic reticulum calcium concentration

  % --- Ca Buffer Concentrations ---
  states(20) = 0.00774891940458; % bc;
  states(21) = 0.01216306138938; % bd;
  states(22) = 12.83478408872024; % bs;
  states(23) = 0.06365082846742; % bsl;

  % --- Membrane potential ---
  states(24) = -79.58430038484185; % V_m;

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