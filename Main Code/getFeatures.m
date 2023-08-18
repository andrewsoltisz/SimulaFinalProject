function [iTrial_flag, feature_outputs_mean, feature_outputs_std, feature_names] = getFeatures(time, state_vals, state_names, stim_startIdx)
    
    % Get Indices of Useful State Variables
    V_idx = matches(state_names, ["V", "V_m"]); % membrane potential (mV) V_m
    Ca_i_idx = matches(state_names, ["ca_i", "cai", "cc"]); % cytosolic [Ca2+] (nM?) cai, cc
    
    percent = [0.3, 0.5, 0.9]; % APDs calculated
    percent_ca = [0.5]; % CaT decay

    %% Feature calculation
    % Names for calculated features
    feature_names = {'Maximum upstroke', 'Minimum diastolic potential', ...
        'AP amplitude', 'APD30%', 'APD50%', 'APD90%', ...
        'Maximum departure velocity', 'Diastolic calcium concentration', ...
        'CaT amplitude', 'CaT time to peak', 'CaT decay 50%', 'CaT tau'};
    % Intialize values and preallocate matrices
    feature_outputs = zeros(length(stim_startIdx),length(feature_names));
    apd = zeros(1, length(percent));
    cat = zeros(1, length(percent_ca));
    
    % Add last index to capture start and end indicies of all beats
    stim_startIdx(end+1) = length(time);

    for iBeat = 1:length(stim_startIdx)-1
        voltage_vals = state_vals(stim_startIdx(iBeat):stim_startIdx(iBeat+1)-1,V_idx);
        time_output = time(stim_startIdx(iBeat):stim_startIdx(iBeat+1)-1);

        %% Voltage feature calculation
        % Maximum upstroke velocity
        [vmax, vmax_idx] = max(voltage_vals);
        [dVdt_max, dVdt_max_idx] = max(diff(voltage_vals(1:vmax_idx,1))./diff(time_output(1:vmax_idx,1)));

        % Minimum diastolic potential
        [vmin, ~] = min(voltage_vals);

        % Action potential amplitude
        ap_amplitude = vmax - vmin;

        % Action potential duration
        for i = 1:length(percent)
            v_apd_end = vmax - ap_amplitude * percent(i);
            v_apd_end_idx = find(voltage_vals(vmax_idx:end) < v_apd_end);
            apd(1,i) = time_output(v_apd_end_idx(1)+vmax_idx, 1) - time_output(dVdt_max_idx, 1);
        end

        %% Calcium feature calculation
        calcium_vals = state_vals(stim_startIdx(iBeat):stim_startIdx(iBeat+1)-1,Ca_i_idx);

        % Maximum departure velocity
        [camax, camax_idx] = max(calcium_vals);
        [dcadt_max, dcadt_max_idx] = max(diff(calcium_vals(1:camax_idx,1))./diff(time_output(1:camax_idx,1)));

        % Diastolic calcium concentration
        [camin, camin_idx] = min(calcium_vals);

        % CaT amplitude
        ca_amplitude = camax - camin;

        % Time to CaT peak
        [~, camin_initial_idx] = min(calcium_vals(1:camax_idx));
        ca_ttp = time_output(camax_idx) - time_output(camin_initial_idx);

        % CaT decay
        for i = 1:length(percent_ca)
            cat_end = camax - ca_amplitude * percent_ca(i);
            cat_end_idx = find(calcium_vals(camax_idx:end) < cat_end);
            cat(1,i) = time_output(cat_end_idx(1)+camax_idx, 1) - time_output(dcadt_max_idx, 1);
        end

        % Time constant of CaT decay (tau)
        cat_end = camax - ca_amplitude * (1-1/exp(1));
        cat_end_idx = find(calcium_vals(camax_idx:end) < cat_end);
        tauFall = time_output(cat_end_idx(1)+camax_idx, 1) - time_output(dcadt_max_idx, 1);

        feature_outputs(iBeat,:) = [dVdt_max, vmin, ap_amplitude, apd, ...
        dcadt_max, camin, ca_amplitude, ca_ttp, cat, tauFall]; 
    end

    %% Calculate average and standard deviation
    feature_outputs_mean = mean(feature_outputs,1);
    feature_outputs_std = std(feature_outputs,0,1);
    
    % check for std variance and set flags
    APD_90_std_thresh = 0.1; 
    if feature_outputs_std(strcmp('APD90%',feature_names))>APD_90_std_thresh*feature_outputs_mean(strcmp('APD90%',feature_names))
        iTrial_flag = false; % this trial is bad
        feature_outputs_mean = zeros(1,length(feature_names));
        feature_outputs_std = zeros(1,length(feature_names));
        feature_names = cell(1,length(feature_names));
    else
        iTrial_flag = true;
    end  
end