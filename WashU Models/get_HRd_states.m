function [states, statenames] = get_HRd_states(varargin)
%initial conditions for state variables

states = zeros(29,1);
statenames = cell(29,1);

temp = load('HRD2004INI.mat');
if isempty(varargin)
    states = temp.X0(4,2:end)'; 
elseif length(varargin)==1
    CL = int(varargin{1});
    states = temp.X0(find(temp.X0(:,1)==CL))';
else 
    disp('varargin to @get_ORd_states must contain only one input')
end
statenames{1} = 'V';
statenames{2} = 'H';
statenames{3} = 'm';
statenames{4} = 'J';
statenames{5} = 'd';
statenames{6} = 'f';
statenames{7} = 'xr';
statenames{8} = 'ca_i'; % intracellular calcium concentration 
statenames{9} = 'na_i'; % intracellular sodium concentration
statenames{10} = 'k_i';
statenames{11} = 'ca_jsr'; % junctional sarcoplasmic reticulum calcium concentration
statenames{12} = 'ca_nsr'; % network sarcoplasmic reticulum calcium concentration
statenames{13} = 'xs';
statenames{14} = 'xs2';
statenames{15} = 'ydv';
statenames{16} = 'ydv2';
statenames{17} = 'zdv';
statenames{18} = 'fca';
statenames{19} = 'fca2';
statenames{20} = 'f2';
statenames{21} = 'dp';
statenames{22} = 'Ctrap';
statenames{23} = 'ro';
statenames{24} = 'ri';
statenames{25} = 'cass';
statenames{26} = 'mL';
statenames{27} = 'hL';
statenames{28} = 'AA';
statenames{29} = 'CL_i';