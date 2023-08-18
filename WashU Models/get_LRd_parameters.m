function [parameters, parameter_names] = get_LRd_parameters()

constantsLRd;
clearvars -except data 

parameters = zeros(81,1);
parameter_names = cell(81,1);
parameters(1)=4.5; parameter_names{1} = 'k_o';
parameters(2)=140; parameter_names{2} = 'na_o';
parameters(3)=1.8; parameter_names{3} = 'ca_o';

names_temp = fieldnames(data);
for i = 1:length(names_temp)
    parameter_names{i+3}=names_temp{i};
    parameters(i+3)=getfield(data,names_temp{i});
end

parameter_names{find(strcmp(parameter_names, 'GK1max'))} = 'g_K1';
parameter_names{find(strcmp(parameter_names, 'GKrmax'))} = 'g_Kr';
parameter_names{find(strcmp(parameter_names, 'GKsmax'))} = 'g_Ks';
parameter_names{find(strcmp(parameter_names, 'pca'))} = 'g_CaL';
parameter_names{find(strcmp(parameter_names, 'GNa'))} = 'g_Na';
parameter_names{find(strcmp(parameter_names, 'gcab'))} = 'g_bCa';
parameter_names{find(strcmp(parameter_names, 'ibarnak'))} = 'g_NaK';
parameter_names{find(strcmp(parameter_names, 'iupbar'))} = 'J_SERCA_bar';
parameter_names{find(strcmp(parameter_names, 'alpha_Rel'))} = 'K_RyR';
parameter_names{end} = 'g_NaCa';

% parameters(18)=50; % changing stimulation start to 50ms
parameters(end)=1; % vNCX parameter introduced for scaling 
