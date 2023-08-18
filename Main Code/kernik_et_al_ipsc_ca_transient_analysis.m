function [amplitude_avg, T2Peak_avg, TransDur_avg, tau_decay_avg, minca_loc, diastolic_ca_avg] = kernik_et_al_ipsc_ca_transient_analysis(Time, Cai)
% Analysis of CaT properties

if length(Cai) > 3000
    [minval_ca] = min(Cai(end-3000:end));
    [maxval_ca] = max(Cai(end-3000:end));
    [~, i] = min(Cai(end-3000:end-2000));
    
elseif length(Cai) > 1000
    [minval_ca] = min(Cai(end-1000:end));
    [maxval_ca] = max(Cai(end-1000:end));
    
else
    [minval_ca] = min(Cai);
    [maxval_ca] = max(Cai);
end
start_position = 105;

check_amp = maxval_ca-minval_ca;

count_peak = 1;
count_min = 1;
peak_find = 0; %Look for a peak start (initial conditions are at rest (MDP))

peak = zeros(1,5);
peak_loc = zeros(1,5);
time_peak = zeros(1,5);
minca = zeros(1,5);
minca_loc = zeros(1,5);
time_minca = zeros(1,5);
amplitude = zeros(1,5);

dCadt = (Cai(2:end)-Cai(1:end-1))./(Time(2:end)-Time(1:end-1));
Time_dCa = Time(2:end);
Ca_acc = (dCadt(2:end)-dCadt(1:end-1))./(Time_dCa(2:end)-Time_dCa(1:end-1));

% figure;
% plot(Time_dCa(2:end), Ca_acc)

for i = start_position:length(Time)-105
    if Cai(i)>=Cai(i-2) && Cai(i)>=Cai(i+2) && peak_find==0 && count_peak<=5  && Cai(i)>maxval_ca-check_amp*.4    %Find peak AP by finding max value after an MDP (goes from positive to negative slope)
        %Time(i+start_position)
        if count_peak <=0
            break
        end
        peak(count_peak)=Cai(i);
        peak_loc(count_peak)=i;
        time_peak(count_peak)=Time(i);
        peak_find=1;
        count_peak=count_peak+1;
    end
    if count_peak>1 && peak_find==1 && count_min<=5
        
        if Cai(i)>peak(count_peak-1)
            peak(count_peak-1)=Cai(i);
            peak_loc(count_peak-1)=i;
            time_peak(count_peak-1)=Time(i);
        end
    end

    if count_peak>1
        if Cai(i)<= Cai(i+1) && (mean(Ca_acc(i:i+3))>=.05*(max(Ca_acc)) || mean(Ca_acc(i:i+3))>=1e3*(mean(Ca_acc))) && peak_find==1 && count_min<=5 && Cai(i)<minval_ca+check_amp*.2 && Time(i)>time_peak(count_peak-1)+200%CaTrans(i)<(max(CaTrans)-.00001) %Find MDP by min value after a peak
            if max(Cai(i:i+100))>=.2*max(Cai)||abs(Time(i)-Time(i+100))<75
                %.05*(max(Ca_acc)) %wprked when anaylsing Wu lab traces
                %      && peak_find==1 && CaTrans(i)<= CaTrans(i-1) &&
                if count_min <=0
                    break
                end
                minca(count_peak-1)=min(Cai(i-100:i+100));
                minca_loc(count_peak-1)=i;
                time_minca(count_peak-1)=Time(i);
                peak_find=0;
                if count_peak>5
                    count_min=count_peak;
                end
            end
        end
        
    end
    
    if count_peak>5 && Cai(i)<maxval_ca-check_amp*.5 % find 5 peaks
        break
    end
    
end

TransDur_ca = zeros(1,count_peak);
TransDur_time=zeros(1,count_peak);
TransDur_loc=zeros(1,count_peak);
TransDur = zeros(1,count_peak);
T2Peak = zeros(1,count_peak);

count = 1;
%count_taucalc=0;
tau_decay_ca = zeros(1,count_peak); tau_decay_time=zeros(1,count_peak); tau_decay_loc=zeros(1,count_peak);
tau_decay = zeros(1,count_peak);

if count_peak == 2
    [minca(1), i] = min(Cai(1:peak_loc(1)));
    %minca(count_peak-1)=CaTrans(i);
    minca_loc(1) = i;
    time_minca(1) = Time(i);
end
if time_peak(1) < time_minca(1)
    time_peak = time_peak(2:end);
    peak = peak(2:end);
    peak_loc = peak_loc(2:end);    
end

for i = start_position:length(Time)-1

    if (count <=length(peak_loc(peak_loc~=0)) && length(minca_loc)>=count+1) %APD measurements must be between a found peak and MDP poitn
        
        amplitude(count) = peak(count)-minca(count); %Amplitude is the voltage change between peak and MDP
        TransDur_end = peak(count)-.9*(amplitude(count)); %APD90
        TransDur_tau = peak(count)-0.632120558828558*(amplitude(count)); %tau decay
        
        if i>peak_loc(count) && (i < minca_loc(count+1)|| minca_loc(count+1)==0) && Cai(i)>TransDur_end   %Find Vm of 90% amplitude between a peak and MDP during repolarization
            TransDur_ca(count) = Cai(i);
            TransDur_time(count) = Time(i);
            TransDur_loc(count) = i;
            
            TransDur(count) = Time(i)-time_minca(count);
            
        end   
        
        if i>peak_loc(count) && (i < minca_loc(count+1)|| minca_loc(count+1)==0 ) && Cai(i)>TransDur_tau   %Find Vm of 90% amplitude between a peak and MDP during repolarization
            tau_decay_ca(count) = Cai(i);
            tau_decay_time(count) = Time(i);
            tau_decay_loc(count) = i;
            
            tau_decay(count) = Time(i)-time_peak(count);
            
        end
        
        if i<max(peak_loc) && T2Peak(count)==0
            T2Peak(count) = time_peak(count)-time_minca(count);
        end
        if i>minca_loc(count+1)&& count_peak>2 %start analyzing next AP, after you pass MDP at end of previous.
            count = count+1;
        end
        
    end
    
    %if i>max(peak_loc) && i>max(minca_loc) %stop analyzing after passing segment analyzed in peak/mdp for loop.
    %    break
    %end
end

%%
amplitude = amplitude(amplitude~=0);
T2Peak = T2Peak(T2Peak~=0);
tau_decay = tau_decay(tau_decay~=0);
TransDur = TransDur(TransDur~=0);
minca = minca(minca~=0);

amplitude_avg = mean(amplitude);
T2Peak_avg = mean(T2Peak);
tau_decay_avg = mean(tau_decay);
TransDur_avg = mean(TransDur);
diastolic_ca_avg = mean(minca);

%% Plot AP with morphology markers
peak = peak(peak_loc~=0);
time_peak = time_peak(peak_loc~=0);
minca = minca(minca_loc~=0);
minca_loc = minca_loc(minca_loc~=0);
time_minca = time_minca(minca_loc~=0);
TransDur_ca = TransDur_ca(TransDur~=0);
TransDur_time = TransDur_time(TransDur~=0);
tau_decay_ca = tau_decay_ca(tau_decay~=0);
tau_decay_time = tau_decay_time(tau_decay~=0);

%% Plot Ca Transient with markers
% figure1=figure;
% axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',18);
% hold(axes1,'all');
% %plot1=plot(Time, CaTrans, 'k', 'Parent',axes1,'LineWidth',5);
% plot1=plot(Time, CaTrans, 'k', 'Parent',axes1,'LineWidth',2);
% hold on
% plot2=plot( time_minca, minca,time_peak, peak, TransDur_time, TransDur_ca, tau_decay_time, tau_decay_ca,'Parent',axes1,'LineWidth',2,'LineStyle','none', 'MarkerSize',30,'Marker','.');
% xlabel('Time (ms)','FontSize',18);
% ylabel('Ca2+ (mM)','FontSize',18);
% hold off

end