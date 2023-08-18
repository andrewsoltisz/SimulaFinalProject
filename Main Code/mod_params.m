function[parameters,parameter_names] = mod_params(parameters,parameter_names,scaling_factors,mod_parameter_names)
    
    SF = scaling_factors;
    mods = mod_parameter_names;
    
    for i =1:length(mods)
        parameters(strcmp(parameter_names,mods{i})) = SF(i)*parameters(strcmp(parameter_names,mods{i}));
    end
    
    
    % mod_parameter_names can in principal be the same as parameter_names,
    % and is in fac the same when we do a Gaussian Convolution to our
    % random matrix.