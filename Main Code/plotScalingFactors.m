function plotScalingFactors(modParam_scaling, modParam_names, nModParams, savePath)    
%  plot scaling factor distributions

    fig = {newFigure('scalingFactors')};
    for iModParam = 1:nModParams
        subplot(4, ceil(nModParams/4), iModParam);
        iScalingDistribution = modParam_scaling(:,iModParam);
        histogram(iScalingDistribution, 'FaceColor', 'blue');
        title(sprintf("%s: Mean = %.3f, Stdev = %.3f", modParam_names{iModParam}, mean(iScalingDistribution), std(iScalingDistribution)), 'Interpreter', 'none');
        grid on;
    end
    sgtitle('Scaling Factor Distributions', 'Interpreter', 'none');
    saveFigures(fig, savePath);

end