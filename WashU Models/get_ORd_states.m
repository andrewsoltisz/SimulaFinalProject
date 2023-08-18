function [states, statenames] = get_ORd_states()
%initial conditions for state variables
states = zeros(41,1);
statenames = cell(41,1);

states(1)=-87; statenames{1} = 'V';
states(2)=7; statenames{2} = 'na_i'; % intracellular sodium concentration
states(3)=145; statenames{3} = 'nass';
states(4)=145; statenames{4} = 'ki';
states(5)=145; statenames{5} = 'kss';
states(6)=1.0e-4; statenames{6} = 'ca_i'; % intracellular calcium concentration
states(7)=1.0e-4; statenames{7} = 'cass'; 
states(8)=1.0; statenames{8} = 'ca_nsr'; % network sarcoplasmic reticulum calcium concentration
states(9)=1.0; statenames{9} = 'ca_jsr'; % junctional sarcoplasmic reticulum calcium concentration
states(10)=0; statenames{10} = 'm';
states(11)=1; statenames{11} = 'hf';
states(12)=1; statenames{12} = 'hs';
states(13)=1; statenames{13} = 'j';
states(14)=1; statenames{14} = 'hsp';
states(15)=1; statenames{15} = 'jp';
states(16)=0; statenames{16} = 'mL';
states(17)=1; statenames{17} = 'hL';
states(18)=1; statenames{18} = 'hLp';
states(19)=0; statenames{19} = 'a';
states(20)=1; statenames{20} = 'iF';
states(21)=1; statenames{21} = 'iS';
states(22)=0; statenames{22} = 'ap';
states(23)=1; statenames{23} = 'iFp';
states(24)=1; statenames{24} = 'iSp';
states(25)=0; statenames{25} = 'd';
states(26)=1; statenames{26} = 'ff';
states(27)=1; statenames{27} = 'fs';
states(28)=1; statenames{28} = 'fcaf';
states(29)=1; statenames{29} = 'fcas';
states(30)=1; statenames{30} = 'jca';
states(31)=0; statenames{31} = 'nca';
states(32)=1; statenames{32} = 'ffp';
states(33)=1; statenames{33} = 'fcafp';
states(34)=0; statenames{34} = 'xrf';
states(35)=0; statenames{35} = 'xrs';
states(36)=0; statenames{36} = 'xs1';
states(37)=0; statenames{37} = 'xs2';
states(38)=1; statenames{38} = 'xk1';
states(39)=0; statenames{39} = 'Jrelnp';
states(40)=0; statenames{40} = 'Jrelp';
states(41)=0; statenames{41} = 'CaMKt';