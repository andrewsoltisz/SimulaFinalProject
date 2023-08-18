function [states, varargout] = base_model_init_states_rabbit()
  % % Default state values for the rabbit base model from:
  % % A Tveito, KH Jaeger, MM Maleckar, WR Giles, and S Wall (2020) 
  % % "Computational translation of drug effects from animal experiments to
  % % human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  % --- Default initial state values --- 
  states = zeros(24, 1);

  % --- I_Na ---
  states(1) = 0.00338797720195; % m;
  states(2) = 0.65360268198698; % j;

  % --- I_NaL ---
  states(3) = 0.00040326523095; % mL;
  states(4) = 0.30802877043842; % hL;

  % --- I_Kr ---
  states(5) = 0.00000744558253; % Xr1;
  states(6) = 0.47043743722699; % Xr2;

  % --- I_Ks ---
  states(7) = 0.00371567254729; % x_Ks;

  % --- I_to ---
  states(8) = 0.90351918425265; % q;
  states(9) = 0.00380750184522; % r;

  % --- I_CaL ---
  states(10) = 0.00000263473214; % d;
  states(11) = 0.99549163805515; % f;
  states(12) = 0.00707626208877; % f_Ca_B;

  % --- I_f ---
  states(13) = 0.25905409322700; % xf;

  % --- RyRs ---
  states(14) = 0.99999999999084; % r_RyR;

  % --- Ca Concentrations ---
  states(15) = 0.12516517727836; % cn;
  states(16) = 0.00003818940023; % cc;
  states(17) = 0.00004774203012; % cd;
  states(18) = 0.00005906286336; % csl;
  states(19) = 0.12489533420468; % cs;

  % --- Ca Buffer Concentrations ---
  states(20) = 0.00342697087250; % bc;
  states(21) = 0.00570376203590; % bd;
  states(22) = 4.35177282725851; % bs;
  states(23) = 0.03412908313068; % bsl;

  % --- Membrane potential ---
  states(24) = -82.00033520583297; % V_m;

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
    state_names{15} = 'cn';
    state_names{16} = 'cc';
    state_names{17} = 'cd';
    state_names{18} = 'csl';
    state_names{19} = 'cs';

    % --- Ca Buffer Concentrations ---
    state_names{20} = 'bc';
    state_names{21} = 'bd';
    state_names{22} = 'bs';
    state_names{23} = 'bsl';

    % --- Membrane potential ---
    state_names{24} = 'V_m';
    varargout(1) = {state_names};
  end
end