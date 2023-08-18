function [states, varargout] = base_model_init_states_zebrafish()
  % % Default state values for the zebrafish base model from:
  % % A Tveito, KH Jaeger, MM Maleckar, WR Giles, and S Wall (2020) 
  % % "Computational translation of drug effects from animal experiments to
  % % human ventricular myocytes", Scientific Reports 10:10537.
  % % doi: 10.1038/s41598-020-66910-0

  % --- Default initial state values --- 
  states = zeros(26, 1);

  % --- I_Na ---
  states(1) = 0.00801076644517; % m;
  states(2) = 0.48988433950339; % j;

  % --- I_NaL ---
  states(3) = 0.00093207234131; % mL;
  states(4) = 0.22022019095546; % hL;

  % --- I_Kr ---
  states(5) = 0.00002789850495; % Xr1;
  states(6) = 0.44961788756812; % Xr2;

  % --- I_Ks ---
  states(7) = 0.00500216122532; % x_Ks;

  % --- I_to ---
  states(8) = 0.87178716972413; % q;
  states(9) = 0.00475528106299; % r;

  % --- I_CaL ---
  states(10) = 0.00000529949807; % d;
  states(11) = 0.99268985625777; % f;
  states(12) = 0.01402239994910; % f_Ca_B;

  % --- I_f ---
  states(13) = 0.25605157362294; % xf;

  % --- RyRs ---
  states(14) = 1.0; % r_RyR;

  % --- Ca Concentrations ---
  states(15) = 0.53784229586268; % cn;
  states(16) = 0.00009460341039; % cc;
  states(17) = 0.00010004578432; % cd;
  states(18) = 0.00015133403529; % csl;
  states(19) = 0.53732106937068; % cs;

  % --- Ca Buffer Concentrations ---
  states(20) = 0.00783433798174; % bc;
  states(21) = 0.01188608421107; % bd;
  states(22) = 12.21882109139196; % bs;
  states(23) = 0.08247049256108; % bsl;

  % --- Membrane potential ---
  states(24) = -77.83279772506836; % V_m;
  
  % --- I_CaT ---
  states(25) = 0.00018448966394; % dCaT;
  states(26) = 0.94854543059327; % dCaT;

  if nargout == 2

    % --- State names --- 
    state_names = cell(26, 1);

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
    
    % --- I_CaT ---
    state_names{25} = 'dCaT';
    state_names{26} = 'dCaT';
    varargout(1) = {state_names};
  end
end