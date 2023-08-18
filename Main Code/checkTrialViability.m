function [iTrial_flag, feature_outputs_mean, feature_outputs_std, feature_names] = checkTrialViability(iTrial_flag, time, state_vals, state_names, param_vals, param_names, stim_startIdx, feature_num)
% This function determines if a a particular parameter trial produces
% viable model behavior over 'lastStims_n' stimulations (beats) applied
% every 'stim_period' ms. The viability flag 'iTrial_flag' starts with a
% state determined by whether the ODE solver was able to converge on a
% solution (true if converged, false otherwise). This function then
% determines if there are any additional reasons to consider this parameter
% trial as non-viable. The 'time' and 'state_vals' inputs are the time and
% state value outputs produces by a set number of stimulations to the
% model, while 'state_names' identifies which state variables belong to
% which columns in 'state_vals' while the 'state_vals' rows correspond to
% different time points defined by the times in each row of 'time'. This
% function will not be called if the ODE solver did not converge for this
% parameter trial, so there is no need to check for that condition here.

    %% Function Outputs Start Empty. Real Values are Set if Trial is Viable.
    feature_outputs_mean = zeros(1,feature_num);
    feature_outputs_std = zeros(1,feature_num);
    feature_names = cell(1,feature_num);

    try % in case of unexpected failure/errors
        %% Get Indices of Useful State Variables
        V_idx = strcmp(state_names, 'V'); % membrane potential (mV)
        Ca_jSR_idx = strcmp(state_names, 'ca_jsr'); % junctional sarcoplasmic reticulum [Ca2+] (nM?)
        Ca_nSR_idx = strcmp(state_names, 'ca_nsr'); % network sarcoplasmic reticulum [Ca2+] (nM?)
        Ca_i_idx = strcmp(state_names, 'ca_i'); % cytosolic [Ca2+] (nM?)
        Na_i_idx = strcmp(state_names, 'na_i'); % cytosolic [Na+] (nM?)
        vjsr_idx = strcmp(param_names, 'vjsr'); % junctional sarcoplasmic reticulum volume (uL?)
        vnsr_idx = strcmp(param_names, 'vnsr'); % network sarcoplasmic reticulum volume (uL?)

        %% Check if Variables Were Found
        V_found = any(V_idx); % membrane potential (mV)
        Ca_jSR_found = any(Ca_jSR_idx); % junctional sarcoplasmic reticulum [Ca2+] (nM?)
        Ca_nSR_found = any(Ca_nSR_idx); % network sarcoplasmic reticulum [Ca2+] (nM?)
        Ca_SR_found = Ca_jSR_found && Ca_nSR_found; % total sarcoplasmic reticulum [Ca2+]
        Ca_i_found = any(Ca_i_idx); % cytosolic [Ca2+] (nM?)
        Na_i_found = any(Na_i_idx); % cytosolic [Na+] (nM?)
        vjsr_found = any(vjsr_idx); % junctional sarcoplasmic reticulum volume (uL?)
        vnsr_found = any(vnsr_idx); % network sarcoplasmic reticulum volume (uL?)

        %% Troubleshooting Plots (Comment-out for real runs)
    %     if V_found
    %         figure;
    %         plot(time, state_vals(:,V_idx));
    %         xlabel("Time (ms)");
    %         ylabel("Membrane Potential (mV)");
    %     end

        %% Compute/Get Filtering Values
        % Membrane Potential Metrics
        if V_found
            V_max = max(state_vals(:, V_idx)); % membrane potential global max (mV)
            V_min = min(state_vals(:, V_idx)); % membrane potential global min (mV)
            V_amp = V_max - V_min; % max membrane potential amplitude over full time course
        end

        % Indices of start and end times for first and last stimulations
        firstStim_startIdx = 1;
        firstStim_endIdx = stim_startIdx(2) - 1;
        firstStim_idx = firstStim_startIdx:firstStim_endIdx;
        lastStim_startIdx = stim_startIdx(end);
        lastStim_endIdx = numel(time);
        lastStim_idx = lastStim_startIdx:lastStim_endIdx;

        % Sarcoplasmic Reticulum [Ca2+] amplitude change
        if Ca_SR_found && vjsr_found && vnsr_found
            % combine junctional and network SR [Ca]
            vjsr_val = param_vals(vjsr_idx); % volume of junctiona SR
            vnsr_val = param_vals(vnsr_idx); % volume of network SR
            Ca_SR_vals = ((state_vals(:, Ca_jSR_idx)*vjsr_val) + (state_vals(:, Ca_nSR_idx)*vnsr_val)) / (vjsr_val + vnsr_val); % volume weighted average [Ca2+]SR
            Ca_SR_amp_firstStim = max(Ca_SR_vals(firstStim_idx)) - min(Ca_SR_vals(firstStim_idx)); % max amplitude of [Ca2+]SR during first stimulation
            Ca_SR_amp_lastStim = max(Ca_SR_vals(lastStim_idx)) - min(Ca_SR_vals(lastStim_idx)); % max amplitude of [Ca2+]SR during last stimulation
            Ca_SR_ampChange = abs(1 - (Ca_SR_amp_firstStim / Ca_SR_amp_lastStim)); % [Ca2+]SR amplitude percent change between first and last stimulation
        end

        % Cytosolic [Ca2+] amplitude change
        if Ca_i_found
            Ca_i_amp_firstStim = max(state_vals(firstStim_idx, Ca_i_idx)) - min(state_vals(firstStim_idx, Ca_i_idx)); % max amplitude of [Ca2+]i during first stimulation
            Ca_i_amp_lastStim = max(state_vals(lastStim_idx, Ca_i_idx)) - min(state_vals(lastStim_idx, Ca_i_idx)); % max amplitude of [Ca2+]i during last stimulation
            Ca_i_ampChange = abs(1 - (Ca_i_amp_firstStim / Ca_i_amp_lastStim)); % [Ca2+]i amplitude percent change between first and last stimulation
        end

        % Cytosolic [Na+] amplitude change (*** change to average Na over whole time range, and add calc. for end diastolic Na (last 50ms) ***)
        if Na_i_found
            Na_i_amp_firstStim = max(state_vals(firstStim_idx, Na_i_idx)) - min(state_vals(firstStim_idx, Na_i_idx)); % max amplitude of [Na+]i during first stimulation
            Na_i_amp_lastStim = max(state_vals(lastStim_idx, Na_i_idx)) - min(state_vals(lastStim_idx, Na_i_idx)); % max amplitude of [Na+]i during last stimulation
            Na_i_ampChange = abs(1 - (Na_i_amp_firstStim / Na_i_amp_lastStim)); % [Na+]i amplitude percent change between first and last stimulation
        end

        %% Check if Any Filtering Values are Outside of Their Viable Range
        if V_found
            % Check if voltage amplitude couldn't exceed 5mV
            if V_amp < 5 %mV
                iTrial_flag = false; % this trial is bad
                return; % no need to check anything else
            end

            % Check if min voltage was greater than -20mV
            if V_min > -20 %mV
                iTrial_flag = false; % this trial is bad
                return; % no need to check anything else
            end
        end

        % Check if peak amplitude change of [Ca2+]SR, [Ca2+]i, or [Na+]i exceeds threshold value
        thresholdChange = 0.02; % percent
        if Ca_i_found
            if Ca_i_ampChange > thresholdChange
                iTrial_flag = false; % this trial is bad
                return; % no need to check anything else
            end
        end
        if Na_i_found
            if Na_i_ampChange > thresholdChange
                iTrial_flag = false; % this trial is bad
                return; % no need to check anything else
            end
        end
        if Ca_SR_found
            if Ca_SR_ampChange > thresholdChange
                iTrial_flag = false; % this trial is bad
                return; % no need to check anything else
            end
        end
    catch
        iTrial_flag = false; % this trial is bad
        return; % no need to check anything else
    end
    
    %% Get Output Features
    if iTrial_flag == true
        try
            [iTrial_flag, feature_outputs_mean, feature_outputs_std, feature_names] = ...
            getFeatures(time, state_vals, state_names, stim_startIdx);
        catch
            iTrial_flag = false; % this trial is bad
            return; % no need to check anything else
        end
    end
    
    % All of the inputs to this function should be sufficient to calculate
    % everything you'll need for feature extraction. For example, parameter
    % values and names are input here so there is no need to generate them
    % again by calling the model's parameter function.
    
    % The input variable 'stim_startIdx' is a vector where each element is
    % the index to the start time for each beat. So, for example,
    % time(stim_startIdx(1)) would return the start time of the first beat,
    % and time(stim_startIdx(2)) would return the start time of the second
    % beat, and so on... Additionally,
    % state_vals(stim_startIdx(1):stim_startIdx(2)-1,V_idx) would return
    % all the voltage values for the first beat, and
    % state_vals(stim_startIdx(2):stim_startIdx(3)-1,V_idx) would return
    % all the voltage values for the second beat, and so on...

end