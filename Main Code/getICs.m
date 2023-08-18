function [modParam_vals, IC_vals, IC_names, trial_flags, feature_outputs_mean, feature_outputs_std, feature_names] = getICs(nStims, stim_freq, modParam_stdev, modParam_scaling, modParam_names, nTrials, model, solveMode)
% This function generates the initial conditions for the model defined by
% the function handle 'modelFun' and its associated parameter and state
% value function handles 'paramFun' and 'stateFun', respectively. Initial
% conditions are defined as the steady state (final) values for all state
% variables following 'nStims' stimulations applied at a frequency of
% 'stim_freq'. The function outputs the matrix 'IC_vals' which contains the
% steady state values for all state variables pertaining to each
% combination of perturbed parameters in 'param_trialVals', the names of
% these variables in 'IC_names', as well as the parameter values
% 'modParam_vals' used to generate the ICs. The IC values in 'IC_vals' for
% each trial are contained in each row and there are as many columns as
% state variables ('nTrials' rows by 'nStates' columns).

    %% Rename Inputs
    model_paramFun = model{1};
    model_stateFun = model{2};
    model_rhsFun   = model{3};
    model_name     = model{4};

    %% Modify Scaling Parameters
    [param_vals_baseline, param_names] = model_paramFun();
    [param_vals_scaled, modParam_baseline, modParam_vals] = modifyParams(param_vals_baseline, param_names, modParam_scaling, modParam_names, nTrials);
%     plotModifiedParams(param_vals_scaled, param_names, modParam_baseline, modParam_names, numel(modParam_names), model_name); 
    
    %% Get Initial Conditions
    IC_file = ['SAVED_', model_name, '_IC.mat'];
    if isfile(IC_file)
        fprintf("%s model initial conditions loaded.\n", model_name);
        load(IC_file);
    else
        % Retrieve Baseline State Values
        [state_vals_baseline, state_names] = model_stateFun();
        nStates = numel(state_vals_baseline);

        % Simulation Parameters
        stepSize = 2e-2; % ms
        options = odeset('RelTol', 1e-06, 'MaxStep', 0.2, 'InitialStep', stepSize);
        stim_period = 1e3 / stim_freq; % ms
        time_range = [0, stim_period];
        
        % Determine Solve-Mode and Determine ICs
        switch solveMode
            case "parallel"
                [IC_vals, runTime] = solveParallel();
            case "plot"
                [IC_vals, runTime] = solveWithPlotting();
            otherwise
                error("Solve mode '%s' is not supported.", solveMode);
        end
        IC_names = state_names';
        fprintf("%.3fs\n", toc(runTime));
        save(IC_file, 'modParam_vals', 'IC_vals', 'IC_names', 'trial_flags', ...
            'feature_outputs_mean', 'feature_outputs_std', 'feature_names', '-v7.3');
    end

    %% Functions to Solve for ICs
    function [IC_vals, runTime] = solveParallel()
    % Solve model using parallel processing
    
        if isempty(gcp('nocreate')) % No parallel pool running 
            parpool; % Startup parallel processing pool
        end
        
        fprintf("%s model: Determining ICs... ", model_name);
        
        % Setup Data-Stream for parfor Waitbar
        D = parallel.pool.DataQueue;
        waitbarHandle = waitbar(0, 'Determining initial conditions...', 'Name', sprintf('%s Model: nTrials = %i, nStims = %i', model_name, nTrials, nStims));
        afterEach(D, @nUpdateWaitbar);
        iIteration = 1;

        % Intialize values and preallocate matrices
        IC_vals = zeros(nTrials, nStates);
        trial_flags = true(nTrials, 1); % true = good trial, false = bad trial
        lastStims_n = 10; % Number of last stimulations to perform trial exlusion filtering on
        lastStims_iFirst = nStims - lastStims_n;
        feature_num = 12; %modify accordingly to feature number
        feature_outputs_mean = zeros(nTrials, feature_num);
        feature_outputs_std = zeros(nTrials, feature_num);
        feature_names = cell(nTrials, feature_num);
        
        % Iteratively solve the model for each parameter trial and store ICs
        runTime = tic;
        parfor iTrial = 1:nTrials  
            state_vals_output = state_vals_baseline';
            time_lastN = 0; % Start time for lastStims_n 
            state_vals_lastN = []; % state variable values for lastStims_n 
            stim_startIdx_lastN = nan(1, lastStims_n); % index for start of each stimulation for last n stims
            for iStim = 1:nStims
                % Solve model for one stimulation
                try
                    [time_output, state_vals_output] = ode15s(model_rhsFun, time_range, state_vals_output(end, :), options, param_vals_scaled(iTrial,:), param_names);
                catch
                    trial_flags(iTrial) = false; % flag this trial (false = bad trial),...
                    break; % skip remaining stimulations, and move onto next trial.
                end
                
                % Check if ODE solver converged, if not, flag this trial and skip to next trial
                if time_output(end) ~= time_range(end) % If the ODE solver was unable to solve out to final time step...
                    trial_flags(iTrial) = false; % flag this trial (false = bad trial),...
                    break; % skip remaining stimulations, and move onto next trial.
                end
                
                % Store last 'lastStims_n' stimulations for trial exlusion filtering
                if iStim > lastStims_iFirst
                    stim_startIdx_lastN(iStim - lastStims_iFirst) = numel(time_lastN);
                    time_lastN = [time_lastN; time_output + time_lastN(end)];
                    state_vals_lastN = [state_vals_lastN; state_vals_output];
                end
            end
            if (nStims >= lastStims_n) && trial_flags(iTrial) % don't check viability if lastN stims weren't collected or if trial has already been exluded
                [trial_flags(iTrial), feature_outputs_mean(iTrial,:), ...
                    feature_outputs_std(iTrial,:), feature_names(iTrial,:)] = ...
                    checkTrialViability(trial_flags(iTrial), time_lastN(2:end), ...
                    state_vals_lastN, state_names, param_vals_scaled(iTrial,:), ...
                    param_names, stim_startIdx_lastN, feature_num);
            end
            IC_vals(iTrial,:) = state_vals_output(end,:); 
            send(D, iTrial); % send waitbar update through parfor data-stream to display
        end
        close(waitbarHandle); 

        % Function for updating parfor waitbar
        function nUpdateWaitbar(~)
            waitbar(iIteration / nTrials, waitbarHandle);
            iIteration = iIteration + 1;
        end
    end

    function [IC_vals, runTime] = solveWithPlotting()
    % Solve model while generating live plots of membrane potential along
    % the way. This solver is not recommended for large values of nTrials.
    
        fprintf("%s model: Determining ICs...", model_name);
        
        % Setup plotting settings
        voltage_idx = strcmp(state_names, 'V');
        legendEntries = strings(1,nTrials);
        figure;
        hold on;
        
        % Intialize values and preallocate matrices
        lastStims_n = 10; % Number of last stimulations to perform trial exlusion filtering on
        lastStims_iFirst = nStims - lastStims_n;
        IC_vals = zeros(nTrials, nStates);
        trial_flags = true(nTrials, 1);  % 1 = good trial, 0 = bad trial
        nIterations = nTrials * nStims;
        iIteration = 0;
        progress = iIteration / nIterations;
        h = waitbar(progress, sprintf('Determining initial conditions: Trial: %i/%i, Stim: %i/%i', 0, nTrials, 0, nStims), 'Name', ['Model: ',model_name]);
        feature_num = 12; %modify accordingly to feature number
        feature_outputs_mean = zeros(nTrials, feature_num);
        feature_outputs_std = zeros(nTrials, feature_num);
        feature_names = cell(nTrials, feature_num);
       
        % Iteratively solve the model for each parameter trial and store ICs
        runTime = tic;
        for iTrial = 1:nTrials 
            time_lastN = 0; % Start time for lastStims_n 
            state_vals_lastN = [];
            stim_startIdx_lastN = nan(1, lastStims_n); % index for start of each stimulation for last n stims
            iTrial_plotColor = [rand, rand, rand];
            legendEntries(iTrial) = "Trial " + string(iTrial);
            time_global = 0;
            state_vals_global = [];
            state_vals_output = state_vals_baseline';
            for iStim = 1:nStims
                iIteration = iIteration + 1;
                progress = iIteration / nIterations;
                waitbar(progress, h, sprintf('Determining initial conditions: Trial: %i/%i, Stim: %i/%i', iTrial, nTrials, iStim, nStims), 'Name', model_name);
                [time_output, state_vals_output] = ode15s(model_rhsFun, time_range, state_vals_output(end, :), options, param_vals_scaled(iTrial,:), param_names);
                
                % Plot membrane potential for each trial
                time_global = [time_global; time_output + time_global(end)];
                state_vals_global = [state_vals_global; state_vals_output];
                voltage_vals = state_vals_global(:, voltage_idx);
                p(iTrial) = plot(time_global(2:end), voltage_vals, 'Color', iTrial_plotColor);
                legend(p,legendEntries(legendEntries ~= ""));
                if iIteration == 1
                    title(sprintf('%s Model: STDEV = %.2f, nTrials = %i, nStims = %i, fStim = %.2fHz', model_name, modParam_stdev, nTrials, nStims, stim_freq));
                    xlabel("Time (ms)");
                    ylabel("Membrane Potential (mV)");
                    grid on;
                end
                drawnow;
                
                % Check if solver converged, if not, flag this trial and skip to next trial
                if time_output(end) ~= time_range(end)
                    trial_flags(iTrial) = false; % flag this trial (false = bad trial)
                    break; % skip remaining stimulations and move onto next trial
                end
                
                % Store last 'lastStims_n' stimulations for trial exlusion filtering
                if iStim > lastStims_iFirst
                    stim_startIdx_lastN(iStim - lastStims_iFirst) = numel(time_lastN);
                    time_lastN = [time_lastN; time_output + time_lastN(end)];
                    state_vals_lastN = [state_vals_lastN; state_vals_output];
                end
            end
            if (nStims >= lastStims_n) && trial_flags(iTrial) % don't check viability if lastN stims weren't collected or if trial has already been exluded
                [trial_flags(iTrial), feature_outputs_mean(iTrial,:), ...
                    feature_outputs_std(iTrial,:), feature_names(iTrial,:)] = ...
                    checkTrialViability(trial_flags(iTrial), time_lastN(2:end), ...
                    state_vals_lastN, state_names, param_vals_scaled(iTrial,:), ...
                    param_names, stim_startIdx_lastN, feature_num);
            end
            IC_vals(iTrial,:) = state_vals_output(end,:);
        end
        close(h);
    end
end