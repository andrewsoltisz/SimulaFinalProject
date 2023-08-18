function [states, statenames] = get_LRd_states(varargin)
%initial conditions for state variables

states = zeros(18,1);
statenames = cell(18,1);

temp = load('Initial_conditions2007.mat');
if isempty(varargin)
    states = temp.X0(35,2:end)'; 
elseif length(varargin)==1
    CL = int(varargin{1});
    states = temp.X0(find(temp.X0(:,1)==CL))';
else 
    disp('varargin to @get_LRd_states must contain only one input')
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
statenames{14} = 'B';
statenames{15} = 'G';
statenames{16} = 'xs2';
statenames{17} = 'Rel';
statenames{18} = 'Over';