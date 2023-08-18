function [feature_baseline_mean, feature_baseline_stdev] = generateBaselineFeatures(model, nStims, stim_freq)

    %% Parse Model Definition 
    model_paramFun = model{1};
    model_stateFun = model{2};
    model_rhsFun   = model{3};
    model_name     = model{4};
    
    %% Check for Previously Generated Values
    featBase_file = ['SAVED_', model_name, '_featBase.mat'];
    if isfile(featBase_file)
        fprintf("%s model feature baseline values loaded.\n", model_name);
        load(featBase_file);
    else
        % Get Baseline/Initial Values for Parameters and State Variables
        [param_vals_baseline, param_names] = model_paramFun();
        [state_vals_baseline, state_names] = model_stateFun();
        
        % Define Simulation Parameters
        stepSize = 2e-2; % ms
        options = odeset('RelTol', 1e-06, 'MaxStep', 0.2, 'InitialStep', stepSize);
        stim_period = 1e3 / stim_freq; % ms
        time_range = [0, stim_period];
        
        % Stimulation and Feature Extraction Settings
        fprintf("%s model: Determining feature baselines... ", model_name);
        lastStims_n = 10; % Number of last stimulations to perform trial exlusion filtering on
        lastStims_iFirst = nStims - lastStims_n;
        feature_num = 12; %modify accordingly to feature number   
        state_vals_output = state_vals_baseline';
        time_lastN = 0; % Start time for lastStims_n 
        state_vals_lastN = []; % state variable values for lastStims_n
        stim_startIdx_lastN = nan(1, lastStims_n); % index for start of each stimulation for last n stims

        % Stimulate the Model
        h = waitbar(0, sprintf('Determining baseline features: Stim: %i/%i', 0, nStims), 'Name', ['Model: ',model_name]);
        runTime = tic;
        for iStim = 1:nStims
            progress = iStim / nStims;
            waitbar(progress, h, sprintf('Determining baseline features: Stim: %i/%i', iStim, nStims), 'Name', ['Model: ',model_name]);
            [time_output, state_vals_output] = ode15s(model_rhsFun, time_range, state_vals_output(end, :), options, param_vals_baseline, param_names);
            if iStim > lastStims_iFirst
                stim_startIdx_lastN(iStim - lastStims_iFirst) = numel(time_lastN);
                time_lastN = [time_lastN; time_output + time_lastN(end)];
                state_vals_lastN = [state_vals_lastN; state_vals_output];
            end
        end
        close(h);
        [~, feature_baseline_mean, feature_baseline_stdev, ~] = checkTrialViability(1, time_lastN(2:end), state_vals_lastN, state_names, param_vals_baseline, param_names, stim_startIdx_lastN, feature_num);
        save(featBase_file, 'feature_baseline_mean', 'feature_baseline_stdev', '-v7.3');
        fprintf("%.3fs\n", toc(runTime));
    end
        
        

end