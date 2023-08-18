function [param_trialVals, modParam_baseline, modParam_vals] = modifyParams(param_vals, param_names, modParam_scaling, modParam_names, nTrials)
% This function creates a population of varied parameter values based on
% the scaling factors in 'modParam_scaling'. The output 'param_trialVals'
% has nTrials rows and a number columns equal to the number of total
% parameters used in the specific model, so the parameter values for a
% specific somulation trial are contained in each row. Only the columns
% pertaining to the modifying parameters get scaled; all other columns
% (parameters) have constant values in each row because they are not being
% scaled/perturbed.

    [~,modParam_idx] = ismember(modParam_names, param_names);
    modParam_baseline = param_vals(modParam_idx)';
    modParam_vals = modParam_baseline .* modParam_scaling;
    param_trialVals = repmat(param_vals', nTrials, 1);
    param_trialVals(:,modParam_idx) = modParam_vals;

end
    
    
    