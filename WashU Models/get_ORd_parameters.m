function [parameters, parameter_names] = get_ORd_parameters()
%initial conditions for state variables
parameters = zeros(31,1);
parameter_names = cell(31,1);

parameters(1)=1; parameter_names{1} = 'celltype';  %endo = 0, epi = 1, M = 2
parameters(2)=8314.0; parameter_names{2} = 'R';
parameters(3)=310.0; parameter_names{3} = 'T';
parameters(4)=96485.0; parameter_names{4} = 'F';
parameters(5)=0.01; parameter_names{5} = 'L';
parameters(6)=0.0011; parameter_names{6} = 'rad';
parameters(7)=1000*pi()*parameters(6)*parameters(6)+parameters(5); parameter_names{7} = 'vcell';
parameters(8)=1000*pi()*parameters(6)*parameters(6)+2*pi()*parameters(6)*parameters(5); parameter_names{8} = 'Ageo';
parameters(9)=parameters(7)*2; parameter_names{9} = 'Acap';
parameters(10)=0.68*parameters(7); parameter_names{10} = 'vmyo';
parameters(11)=0.0552*parameters(7); parameter_names{11} = 'vnsr'; % network sarcoplasmic reticulum volume
parameters(12)=0.0048*parameters(7); parameter_names{12} = 'vjsr'; % junctional sarcoplasmic reticulum volume
parameters(13)=0.02*parameters(7); parameter_names{13} = 'vss';
parameters(14)=0.1908; parameter_names{14} = 'g_K1';
parameters(15)=0.046; parameter_names{15} = 'g_Kr';
parameters(16)=0.0034; parameter_names{16} = 'g_Ks';
parameters(17)=0.02; parameter_names{17} = 'g_to';
parameters(18)=0.00005; parameter_names{18} = 'g_CaL';
parameters(19)=2.5e-8; parameter_names{19} = 'g_bCa';
parameters(20)=75; parameter_names{20} = 'g_Na';
parameters(21)=0.0075; parameter_names{21} = 'g_NaL';
parameters(22)=30; parameter_names{22} = 'g_NaK';
parameters(23)=0.0008; parameter_names{23} = 'g_NaCa';
parameters(24)=0.0005; parameter_names{24} = 'g_pCa';
parameters(25)=0.004375; parameter_names{25} = 'J_SERCA_bar';
parameters(26)=1; parameter_names{26} = 'K_RyR';
parameters(27)=140.0; parameter_names{27} = 'Nao';
parameters(28)=1.8; parameter_names{28} = 'Cao';
parameters(29)=5.4; parameter_names{29} = 'Ko';
parameters(30)=1; parameter_names{30} = 'ODE';
parameters(31)=0; parameter_names{31} = 'Vclamp';