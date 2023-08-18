function [states, varargout] = base_model_init_states_dog()
  % % Default state values for the canine base model from:
  % % A Tveito, KH Jaeger, MM Maleckar, WR Giles, and S Wall (2020) 
  % % "Computational translation of drug effects from animal experiments to
  % % human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  % --- Default initial state values --- 
  states = zeros(24, 1);

  % --- I_Na ---
  states(1) = 0.00101330247289; % m;
  states(2) = 0.81832895816850; % j;

  % --- I_NaL ---
  states(3) = 0.00012984267222; % mL;
  states(4) = 0.48336584886284; % hL;

  % --- I_Kr ---
  states(5) = 0.00000114164325; % Xr1;
  states(6) = 0.49873718562668; % Xr2;

  % --- I_Ks ---
  states(7) = 0.00248189604333; % x_Ks;

  % --- I_to ---
  states(8) = 0.93540710727203; % q;
  states(9) = 0.00281717928533; % r;

  % --- I_CaL ---
  states(10) = 0.00000102457721; % d;
  states(11) = 0.99777059597759; % f;
  states(12) = 0.00668797146277; % f_Ca_B;

  % --- I_f ---
  states(13) = 0.34750343894071; % xf;

  % --- RyRs ---
  states(14) = 0.99999999997792; % r_RyR;

  % --- Ca Concentrations ---
  states(15) = 0.16870006802854; % cn, network sarcoplasmic reticulum calcium concentration
  states(16) = 0.00004311732466; % cc; % cytosolic calcium concentration
  states(17) = 0.00004511348609; % cd;
  states(18) = 0.00006825503428; % csl;
  states(19) = 0.16840409968237; % cs, junctional sarcoplasmic reticulum calcium concentration

  % --- Ca Buffer Concentrations ---
  states(20) = 0.00383339141962; % bc;
  states(21) = 0.00539092855004; % bd;
  states(22) = 5.55581916181986; % bs;
  states(23) = 0.03919823900131; % bsl;

  % --- Membrane potential ---
  states(24) = -87.70510492702755; % V_m;

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