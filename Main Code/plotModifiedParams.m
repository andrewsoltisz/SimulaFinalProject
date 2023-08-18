function plotModifiedParams(param_trialVals, param_names, modParam_baseline, modParam_names, nModParams, model_name)    
% plot modified parameter distributions

    [~,modParam_idx] = ismember(modParam_names, param_names);
    figure;
    for iModParam = 1:nModParams
        subplot(4, ceil(nModParams/4), iModParam);
        iParamDistribution = param_trialVals(:,modParam_idx(iModParam));
        histogram(iParamDistribution);
        title(sprintf("%s: Base = %f, Mean = %f, Stdev = %f", modParam_names{iModParam}, modParam_baseline(iModParam), mean(iParamDistribution), std(iParamDistribution)), 'Interpreter', 'none');
        grid on;
    end
    sgtitle([model_name, ' Model: Trial Parameter Distributions'], 'Interpreter', 'none')

end