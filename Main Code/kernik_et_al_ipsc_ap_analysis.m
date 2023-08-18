function [AP_measures, flag] = kernik_et_al_ipsc_ap_analysis(Time, values, sim_inputs, SA_par)
% Analysis of AP and CaT properties

if length(Time)>10000
    values=values(end-10000:end,:);
    Time=Time(end-10000:end);
end

num_measures=14;
flag=zeros(1,7);
Vm=values(:,1);
options = odeset('MaxStep',1,'InitialStep',2e-5);
if sim_inputs(1)==1
    flag(1)=1;
end

if length(Vm)>1000
    if sim_inputs(1)==0 && max(Vm(500:end))-min(Vm(500:end))< 5 && min(Vm(500:end)) <-20 %Non beating, and rest (less than 5mV change over time course), and rest is less than -20mV
        flag(1)=1;
        sim_inputs(1)=1; %Implement stim
        % [ Y_IC ] = Find_SS( values(end,:), options,model_parameter_inputs );
        [Time1, values1] = ode15s(@kernik_et_al_ipsc_function,[0, 59990],values(end,:), options, sim_inputs, SA_par); %Run with beating for 59 beats, pull ICs after last beat
        [Time, values] = ode15s(@kernik_et_al_ipsc_function,[0, 10e3],values1(end,:), options, sim_inputs, SA_par); %run for analysis
        Vm=values(:,1);
    end
end

    Ca_i=values(:,3);
    Ca_SR=values(:,2);
    Na_i=values(:,4);

if max(Time)< 1 %if run fails
    AP_measures=zeros(1,num_measures);
    flag(5)=1;
elseif length(Ca_i)<6000 || length(Na_i)<6000 %if run less than 4000 time steps
    AP_measures=zeros(1,num_measures);
    flag(5)=1;
elseif max(Vm(500:end))-min(Vm(500:end)) < 5 && flag(1)==0%Non-functional, non-beating cells wiht rest voltage >-20mV
    flag(2)=1;
    AP_measures=zeros(1,num_measures);
else
    Ca_i_min=min(Ca_i(end-3000:end));
    Na_i_min=min(Na_i(end-3000:end));
    Ca_i_min_start=min(Ca_i(1:3000));
    Na_i_min_start=min(Na_i(1:3000));
    Ca_i_max=max(Ca_i(end-3000:end));
    Na_i_max=max(Na_i(end-3000:end));
    Ca_i_max_start=max(Ca_i(1:3000));
    Na_i_max_start=max(Na_i(1:3000));
    Ca_SR_min=min(Ca_SR(end-3000:end));
    Ca_SR_min_start=min(Ca_SR(1:3000));
    Ca_SR_max=max(Ca_SR(end-3000:end));
    Ca_SR_max_start=max(Ca_SR(1:3000));
    
    if flag(1)~=1 && abs(1-(Ca_i_max_start-Ca_i_min_start)/(Ca_i_max-Ca_i_min))>.02 || abs(1-(Ca_SR_max_start-Ca_SR_min_start)/(Ca_SR_max-Ca_SR_min))>.02% || abs(1-(Na_i_max_start-Na_i_min_start)/(Na_i_max-Na_i_min))>.02
        %if not at SS during run (amplitude of Cai, Nai, or Ca_SR changes
        %more than 2% over run
        
        flag(4)=1;
        
    end
  
%     if abs(1-(Ca_i_max_start-Ca_i_min_start)/(Ca_i_max-Ca_i_min))>.02
%         cai_flag=1
%     end
%     
%     if abs(1-(Ca_SR_max_start-Ca_SR_min_start)/(Ca_SR_max-Ca_SR_min))>.02
%         abs(1-(Ca_SR_max_start-Ca_SR_min_start)/(Ca_SR_max-Ca_SR_min))
%         caSR_flag=1
%     end
%     
%     if abs(1-(Na_i_max_start-Na_i_min_start)/(Na_i_max-Na_i_min))>.02
%         na_flag=1
%     end
    
    %Initialize variables:
    count_peak=1;
    count_min=1;
    peak_find=0; %Look for a peak start (initial conditions are at rest (MDP)
    peak=zeros(1,5); peak_loc=zeros(1,5); time_peak=zeros(1,5);
    dv_dt_max=zeros(1,5); dvmax_loc=zeros(1,5);
    time_dvmax=zeros(1,5); vm_dvmax=zeros(1,5);
    mdp=zeros(1,5); mdp_loc=zeros(1,5); time_mdp=zeros(1,5);
    amplitude=zeros(1,5);
    APD90_vm=zeros(1,5); APD90_time=zeros(1,5); APD90_loc=zeros(1,5);
    APD30=zeros(1,5);
    APD40=zeros(1,5);
    APD50=zeros(1,5);
    APD70=zeros(1,5);
    APD80=zeros(1,5);
    APD90=zeros(1,5);
    
    for i=2:length(Time)-1
        if Vm(i)>Vm(i-1) && Vm(i)>Vm(i+1) && peak_find==0 && count_peak<=5  %Find peak AP by finding max value after an MDP (goes from positive to negative slope)
            if Vm(i)>-20 %"real" peaks should have voltage >=-20mV
                if count_peak <=0 %if oscillating is happening, and causes counter to get below zero, set flag for ossilating.
                    flag(6)=1;
                    break
                end
                peak(count_peak)=Vm(i);
                peak_loc(count_peak)=i;
                time_peak(count_peak)=Time(i);
                peak_find=1;
                count_peak=count_peak+1;
                
                if count_min >1 && (count_peak-2)>0
                    if Time(i)-100 < time_peak(count_peak-2) && Vm(i)+20 > mdp(count_min-1) %if peak is foung within 20mV and 100ms of previous peak
                        %reset  min to previous min (ignore min between
                        %true peak and notch
                        count_min=count_min-1;
                        count_peak=count_peak-1;
                    end
                end
                
            elseif sum(peak)>0 %if peak is less than -20mV, ignore previous minima, find new min (hyperpolarization of AP)
                count_min=count_min-1;
                peak_find=1;
            end
            
        end
        
        if Vm(i)< Vm(i-1) && Vm(i)< Vm(i+1) && peak_find==1 && count_min<=5  %Find MDP by min value after a peak
            if count_min <=0 %if oscillating is happening, and causes counter to get below zero, set flag for ossilating.
                flag(6)=1;
                break
            end
            mdp(count_min)=Vm(i);
            mdp_loc(count_min)=i;
            time_mdp(count_min)=Time(i);
            peak_find=0;
            count_min=count_min+1;
            
        end
        if count_peak>5 %find 5 peaks
            break
        end
        
    end
    
    %Define dv/dt for full trace:
    Vm_nminus1=Vm(2:end);
    Vm_n=Vm(1:end-1);
    t_nminus1=Time(2:end);
    t_n=Time(1:end-1);
    v_prime=(Vm_n-Vm_nminus1)./(t_n-t_nminus1);
    v_prime=[v_prime; zeros(length(Time)-length(v_prime),1)];
    %figure
    %plot(Time, v_prime, Time, Vm./(max(Vm)))
    
    count=1;
    
    for i=2:length(Time)-1
        
        if count>1
            if i<peak_loc(count) && i > mdp_loc(count-1) && v_prime(i)>dv_dt_max(count-1)
                %max dV/dT during depolarization, between MDP and peak
                dv_dt_max(count-1)=v_prime(i);
                vm_dvmax(count-1)=Vm(i);
                time_dvmax(count-1)=Time(i);
                dvmax_loc(count-1)=i;
            end
        end
        
        if count <=length(peak_loc(peak_loc~=0)) && count <=length(mdp_loc(mdp_loc~=0)) %APD measurements must be between a found peak and MDP poitn
            amplitude(count)=peak(count)-mdp(count); %Amplitude is the voltage change between peak and MDP
            APD_end=peak(count)-.9*(amplitude(count)); %APD90
            APD_end30=peak(count)-.3*(amplitude(count));
            APD_end40=peak(count)-.4*(amplitude(count));
            APD_end50=peak(count)-.5*(amplitude(count));
            APD_end70=peak(count)-.7*(amplitude(count));
            APD_end80=peak(count)-.8*(amplitude(count));
            
            
            if i>peak_loc(count) && i < mdp_loc(count) && Vm(i)>APD_end  %Find Vm of 90% amplitude between a peak and MDP during repolarization
                APD90_vm(count)=Vm(i);
                APD90_time(count)=Time(i);
                APD90_loc(count)=i;
                APD90(count)=Time(i);
            end
            if i>peak_loc(count) && i < mdp_loc(count) && Vm(i)>APD_end30   %Find Vm of 30% amplitude between a peak and MDP during repolarization
                APD30(count)=Time(i);
            end
            if i>peak_loc(count) && i < mdp_loc(count) && Vm(i)>APD_end40  %Find Vm of 40% amplitude between a peak and MDP during repolarization
                APD40(count)=Time(i);
            end
            if i>peak_loc(count) && i < mdp_loc(count) && Vm(i)>APD_end50    %Find Vm of 50% amplitude between a peak and MDP during repolarization
                APD50(count)=Time(i);
            end
            if i>peak_loc(count) && i < mdp_loc(count) && Vm(i)>APD_end70   %Find Vm of 70% amplitude between a peak and MDP during repolarization
                APD70(count)=Time(i);
            end
            if i>peak_loc(count) && i < mdp_loc(count) && Vm(i)>APD_end80   %Find Vm of 80% amplitude between a peak and MDP during repolarization
                APD80(count)=Time(i);
            end
            
            
            if i>mdp_loc(count) %start analyzing next AP, after you pass MDP at end of previous.
                count=count+1;
            end
        end
        
        if i>max(peak_loc) && i>max(mdp_loc) %stop analyzing after passing segment analyzed in peak/mdp for loop.
            break
        end
    end
    %APD is
    size_APDs=length(dv_dt_max(dv_dt_max~=0));
    if length(APD90(APD90~=0))>size_APDs
        %APD after first peak, doesn't have corresponding dv_dt_max
        APD30=APD30(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD40=APD40(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD50=APD50(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD70=APD70(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD80=APD80(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD90=APD90(2: size_APDs+1)-time_dvmax(1: size_APDs);
    else
        size_APDs=length(APD90(APD90~=0))-1;
        APD30=APD30(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD40=APD40(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD50=APD50(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD70=APD70(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD80=APD80(2: size_APDs+1)-time_dvmax(1: size_APDs);
        APD90=APD90(2: size_APDs+1)-time_dvmax(1: size_APDs);
        
    end
    
    plateau=(APD30-APD40)./(APD70-APD80); %From Ma et al. ventricular like = 2.5 +/- .2, atrial like is 1.1 +/- .1
    
    %dv_dt_max=dv_dt_max(1:length(dv_dt_max(dv_dt_max~=0)));
    
    cycle_length=time_mdp(2)-time_mdp(1);
    BPM=60e3/cycle_length;
    
    if flag(1)==1 && abs(BPM-60)>5 %if beating leads to finding new steady-state with spontaneous beating
        flag(1)=0;
        flag(7)=1;
    end
        
    if length(mdp(mdp~=0))<1 || length(amplitude(amplitude~=0))<1 || length(APD90(APD90~=0))<1 || length(dv_dt_max(dv_dt_max~=0))<1|| length(APD50(APD50~=0))<1|| length(plateau(plateau~=0))<1 %throw it out if something failed.
        flag(2)=1;
        flag(3:6)=zeros(1,4);
        AP_measures=zeros(1,num_measures);
    elseif sum(flag(2:6))~=0
        AP_measures=zeros(1,num_measures);
        %AP_measures=[-1*mdp(1),amplitude(1), APD90(1), dv_dt_max(1), BPM, Ca_i_min, Na_i_min, APD50(1), plateau(1), cycle_length, 0, 0, 0, 0];
        
    else
        [Ca_amplitude, Ca_T2Peak, Ca_TransDur, Ca_tau_decay, minca_loc, diastolic_ca_avg] = kernik_et_al_ipsc_ca_transient_analysis(Time, Ca_i);
        AP_measures = [-1*mdp(1), amplitude(1), APD90(1), dv_dt_max(1), BPM, Ca_i_min, Na_i_min, APD50(1), plateau(1), cycle_length, Ca_amplitude, Ca_T2Peak, Ca_TransDur, Ca_tau_decay];
        
        %Check for alternans-type behaviors
        if std(amplitude(1:size_APDs)) > 1
            flag(3)=1;
        end
    end
    
    %figure
    %plot(Time, Vm,'k', time_dvmax, vm_dvmax,'.', time_peak, peak, '.',time_mdp, mdp, '.', APD90_time, APD90_vm, '.')
end
end