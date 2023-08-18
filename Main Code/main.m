% SSCP 2021, Project 5
% Advisor: Dr. Andrew G. Edwards
% Students: John Dawson, Anna Gams, Ivan Vikram Rajen, Andrew Soltisz
% Project Title: Computational prediction of cardiac electropharmacology - how much does the model matter?

%% Prepare MATLAB Environment

close all; clear; clc; 
totalRunTime = tic;
currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

%% Quick-Access Analysis Settings

nTrials = 1000;          % number of different parameter combinations (target = 1000)
modParam_stdev = 0.3; % (norm) deviation applied to specified trial parameters (target = TBD)
nStims = 500;         % number of stimulations applied to models before ICs are calculated (target = 500)
stim_freq = 1;        % (Hz) frequency of model stimulation (target = 1Hz)
modParam_names = {'g_K1', 'g_Kr', 'g_Ks', 'g_CaL', 'g_bCa', 'g_Na',... % trial parameters to be modified/perturbed
                  'g_NaK', 'g_NaCa', 'J_SERCA_bar', 'K_RyR'};
              
%% Define Models to Analyze
% model = [model, {{parameter_function_handle, state_function_handle, rhs_function_handle, 'modelName'}}]; % example format (comment-out when running code)

% Individual models (code lines) can be commented-out without affecting how the code runs/works
model = cell(0);
model = [model, {{@get_ORd_parameters, @get_ORd_states, @ORd, 'WU Human'}}];
% model = [model, {{@get_LRd_parameters, @get_LRd_states, @LRd, 'WU Guinea Pig'}}];
model = [model, {{@get_HRd_parameters, @get_HRd_states, @HRd, 'WU Canine'}}];
model = [model, {{@base_model_init_parameters_human, @base_model_init_states_human, @base_model_rhs, 'Simula Human'}}];
model = [model, {{@base_model_init_parameters_guinea_pig, @base_model_init_states_guinea_pig, @base_model_rhs, 'Simula Guinea Pig'}}];
model = [model, {{@base_model_init_parameters_dog, @base_model_init_states_dog, @base_model_rhs, 'Simula Canine'}}];

nModels = numel(model);

%% Setup Model Parameters to be Perturbed

scalingFactorFile = 'SAVED_scalingFactors.mat';
if isfile(scalingFactorFile)
    fprintf("Scaling factors loaded.\n");
    load(scalingFactorFile);
else
    fprintf("Scaling factors generated.\n");
    nModParams = numel(modParam_names);
    modParam_scaling = getScalingFactors(modParam_stdev, nModParams, nTrials);
    save(scalingFactorFile, 'modParam_scaling', 'modParam_names', 'nModParams', 'modParam_stdev', 'nTrials', 'nStims', 'stim_freq', '-v7.3');
end

figureFolder_general = [pwd, filesep, 'Figures'];
if exist(figureFolder_general, 'dir') ~= 7
    mkdir(figureFolder_general);
end
plotScalingFactors(modParam_scaling, modParam_names, nModParams, figureFolder_general);

%% Determine Model Initial Conditions (ICs) and Output Features

% solveMode options are "parallel" (for parallel processing) or "plot" (to plot voltage state vs time, NOT RECOMMENDED FOR LARGE nTrials!)
solveMode = "parallel"; 
% solveMode = "plot"; 

modParam_vals = cell(nModels, 1);
IC_vals = cell(nModels, 1);
IC_names = cell(nModels, 1);
trial_flags = cell(nModels, 1);
feature_outputs_mean = cell(nModels, 1);
feature_outputs_std = cell(nModels, 1);
feature_names = cell(nModels, 1);
for iModel = 1:nModels
    [modParam_vals{iModel}, IC_vals{iModel}, IC_names{iModel}, trial_flags{iModel}, ...
        feature_outputs_mean{iModel}, feature_outputs_std{iModel}, feature_names{iModel}] = ...
        getICs(nStims, stim_freq, modParam_stdev, modParam_scaling, modParam_names, nTrials, ... 
        model{iModel}, solveMode);
end

    %% Perform Intra-model Input-Output Sensitivity Analysis and Plot Results
    
    wb_sensativity = waitbar(0, sprintf('Sensitivity analysis (%i/%i): %s',0,nModels,model{1}{4}), 'Name', 'Sensitivity Analysis');
    sensitivity_forward = cell(nModels, 1); 
    sensitivity_reverse = cell(nModels, 1);
    R2_sens_forward = cell(nModels, 1);
    R2_sens_reverse = cell(nModels, 1);
    
    for iModel = 1:nModels
        waitbar(iModel/nModels, wb_sensativity, sprintf('Sensitivity analysis (%i/%i): %s',iModel,nModels,model{iModel}{4}), 'Name', 'Sensitivity Analysis');
        % flagged trials (false/bad) are ommitted from the PLS function input
        X_raw = modParam_vals{iModel}(trial_flags{iModel}, :); 
        Y_raw = feature_outputs_mean{iModel}(trial_flags{iModel}, :);
        X_log = log(X_raw);
        Y_log = log(Y_raw);
        sensitivity_forward{iModel} = PLS_nipals(X_log, Y_log); 
        sensitivity_reverse{iModel} = PLS_nipals(Y_log, X_log);
        R2_sens_forward{iModel} = plotPLS(sensitivity_forward{iModel}, X_raw, Y_raw, Y_log, modParam_names, feature_names{iModel}(1,:), model{iModel}{4}, model{iModel}{4}, 'Sensitivity', 'Forward', 'Raw', figureFolder_general);
        R2_sens_reverse{iModel} = plotPLS(sensitivity_reverse{iModel}, Y_raw, X_raw, X_log, feature_names{iModel}(1,:), modParam_names, model{iModel}{4}, model{iModel}{4}, 'Sensitivity', 'Reverse', 'Raw', figureFolder_general);
    end
    close(wb_sensativity);

    %% Get Model Feature Baseline Values and Normalize Output Features

    feature_baseline_mean = cell(nModels, 1);
    feature_baseline_stdev = cell(nModels, 1);
    feature_outputs_norm = cell(nModels, 1);

    for iModel = 1:nModels
        [feature_baseline_mean{iModel}, feature_baseline_stdev{iModel}] = generateBaselineFeatures(model{iModel}, nStims, stim_freq);
        feature_outputs_norm{iModel} = (feature_outputs_mean{iModel} - feature_baseline_mean{iModel}) ./ feature_baseline_mean{iModel};
    end

    %% Perform Inter-model Translation Analysis and Plot Results

    model_pairs = nchoosek(1:nModels, 2); % determine model pairs
    nPairs = size(model_pairs, 1); % count total number of model pairs
    translation_forward_raw = cell(nPairs, 1); 
    translation_reverse_raw = cell(nPairs, 1);
    translation_forward_norm = cell(nPairs, 1); 
    translation_reverse_norm = cell(nPairs, 1);
    R2_trans_forward_raw = cell(nPairs, 1);
    R2_trans_reverse_raw = cell(nPairs, 1);
    R2_trans_forward_norm = cell(nPairs, 1);
    R2_trans_reverse_norm = cell(nPairs, 1);

    wb_translation = waitbar(0, sprintf('Translation analysis (%i/%i): %s-%s',0,nPairs,model{model_pairs(1,1)}{4},model{model_pairs(1,2)}{4}), 'Name', 'Translation Analysis');
    for iPair = 1:nPairs
        model_1 = model_pairs(iPair, 1); % get the number of the first model (for simpler syntax)
        model_2 = model_pairs(iPair, 2); % get the number of the second model
        waitbar(iPair/nPairs, wb_translation, sprintf('Translation analysis (%i/%i): %s-%s',iPair,nPairs,model{model_1}{4},model{model_2}{4}), 'Name', 'Translation Analysis');
        trial_flags_pair = trial_flags{model_1} & trial_flags{model_2}; % only trials that are not flagged (true/good) in both models are passed to the PLS function

        % Analysis of Raw Values
        X_raw = feature_outputs_mean{model_1}(trial_flags_pair, :);
        Y_raw = feature_outputs_mean{model_2}(trial_flags_pair, :);
        X_log = log(X_raw);
        Y_log = log(Y_raw);
        translation_forward_raw{iPair} = PLS_nipals(X_log, Y_log);
        translation_reverse_raw{iPair} = PLS_nipals(Y_log, X_log);
        R2_trans_forward_raw{iPair} = plotPLS(translation_forward_raw{iPair}, X_raw, Y_raw, Y_log, feature_names{model_1}(1,:), feature_names{model_2}(1,:), model{model_1}{4}, model{model_2}{4}, 'Translation', 'Forward', 'Raw', figureFolder_general);
        R2_trans_reverse_raw{iPair} = plotPLS(translation_reverse_raw{iPair}, Y_raw, X_raw, X_log, feature_names{model_2}(1,:), feature_names{model_1}(1,:), model{model_2}{4}, model{model_1}{4}, 'Translation', 'Reverse', 'Raw', figureFolder_general);

        % Analysis of Normalized Values (BROKEN)
    %     X_norm = feature_outputs_norm{model_1}(trial_flags_pair, :);
    %     Y_norm = feature_outputs_norm{model_2}(trial_flags_pair, :);
    %     X_logNorm = log(X_norm);
    %     Y_logNorm = log(Y_norm);
    %     translation_forward_norm{iPair} = PLS_nipals(X_logNorm, Y_logNorm);
    %     translation_reverse_norm{iPair} = PLS_nipals(Y_logNorm, X_logNorm);
    %     R2_trans_forward_norm{iPair} = plotPLS(translation_forward_norm{iPair}, X_norm, Y_norm, Y_logNorm, feature_names{model_1}(1,:), feature_names{model_2}(1,:), model{model_1}{4}, model{model_2}{4}, 'Translation', 'Forward', 'Normalized', figureFolder_general);
    %     R2_trans_reverse_norm{iPair} = plotPLS(translation_reverse_norm{iPair}, Y_norm, X_norm, X_logNorm, feature_names{model_2}(1,:), feature_names{model_1}(1,:), model{model_2}{4}, model{model_1}{4}, 'Translation', 'Reverse', 'Normalized', figureFolder_general);
    end
    close(wb_translation);

model_names = cellfun(@(m) m{4}, model, 'UniformOutput', false);
plotSummary(R2_sens_forward, R2_sens_reverse,...
            R2_trans_forward_raw, R2_trans_reverse_raw,...
            model_names, model_pairs, figureFolder_general);

%% Cross-validation
% 
% translation_forward_train = cell(nPairs, 1); 
% translation_reverse_train = cell(nPairs, 1);
% translation_forward_test = cell(nPairs, 1); 
% translation_reverse_test = cell(nPairs, 1);
% for iPair = 1:nPairs
%     model_1 = model_pairs(iPair, 1); % get the number of the first model (for simpler syntax)
%     model_2 = model_pairs(iPair, 2); % get the number of the second model
%     % only trials that are not flagged (true/good) in both models are passed to the PLS function
%     trial_flags_pair = trial_flags{model_1} & trial_flags{model_2}; 
%     X_norm = feature_outputs_mean{model_1}(trial_flags_pair, :);
%     Y_norm = feature_outputs_mean{model_2}(trial_flags_pair, :);
%     
%     idx = randperm(length(X_norm)); % generate random indicies
%     [~, columns] = size(idx);
%     lastRow = floor(0.8 *columns); 
%     train_idx = idx(1, 1:lastRow);
%     test_idx = idx(1, lastRow+1:end);
%     
%     X_train = X_norm(train_idx, :); % split into train and test sets
%     Y_train = Y_norm(train_idx, :);
%     
%     X_test = X_norm(test_idx, :);
%     Y_test = Y_norm(test_idx, :);
%     
%     % training
%     translation_forward_train{iPair} = PLS_nipals(X_train, Y_train);
%     translation_reverse_train{iPair} = PLS_nipals(Y_train, X_train);
%     plotPLS(translation_forward_train{iPair}, X_train, Y_train, feature_names{model_1}, feature_names{model_2}(1,:), model{model_1}{4}, model{model_2}{4}, 'translation', 'forward', 'raw');
%     plotPLS(translation_reverse_train{iPair}, Y_train, X_train, feature_names{model_2}, feature_names{model_1}(1,:), model{model_1}{4}, model{model_2}{4}, 'translation', 'reverse', 'raw');
%     
%     % testing %not complete!
%     translation_forward_test{iPair} = PLS_nipals_cross_validation(X_test, Y_test, translation_forward_train{iPair});
%     translation_reverse_test{iPair} = PLS_nipals_cross_validation(X_test, Y_test, translation_reverse_train{iPair});
%     % plots
%     
% end

%% Conclude Analysis

fprintf("\nAnalysis Complete!\nTotal run time: %.3fs\n", toc(totalRunTime));