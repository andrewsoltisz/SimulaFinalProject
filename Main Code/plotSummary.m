function plotSummary(R2_sens_forward, R2_sens_reverse,...
            R2_trans_forward, R2_trans_reverse,...
            model_names, model_pairs, savePath)
        
    fprintf("Generating summary plots... ");
    runTime = tic;
    figures = {};
    
    %% Sensitivity R^2 Values
    
    % Forward
    figures = [figures, {newFigure('R2 Sensitivity')}];
    subplot(2, 1, 1);
    R2_avg = cellfun(@(R2) R2{1}, R2_sens_forward);
    R2_std = cellfun(@(R2) R2{2}, R2_sens_forward);
    hold on;
    b = bar(R2_avg, 'blue');
    errorbar(R2_avg, R2_std, 'k', 'linestyle', 'none');
    set(gca,'XTick',b.XData,'XTickLabel',model_names);
    title('Forward Sensitivity');
    ylabel('R^2');
    ylim([0,1]);
    yticks(0:0.2:1);
    grid on;
    hold off;
    
    % Reverse
    subplot(2, 1, 2);
    R2_avg = cellfun(@(R2) R2{1}, R2_sens_reverse);
    R2_std = cellfun(@(R2) R2{2}, R2_sens_reverse);
    hold on;
    b = bar(R2_avg, 'blue');
    errorbar(R2_avg, R2_std, 'k', 'linestyle', 'none');
    set(gca,'XTick',b.XData,'XTickLabel',model_names);
    title('Reverse Sensitivity');
    ylabel('R^2');
    ylim([0,1]);
    yticks(0:0.2:1);
    grid on;
    hold off;
    
     %% Translation R^2 Values (FIX ME)
     
     pair_names = string(model_names(model_pairs));
     forward_names = pair_names(:,1) + "\rightarrow" + pair_names(:,2);
     reverse_names = pair_names(:,1) + "\leftarrow" + pair_names(:,2);
    
    % Forward
    figures = [figures, {newFigure('R2 Translations')}];
    subplot(2, 1, 1);
    R2_avg = cellfun(@(R2) R2{1}, R2_trans_forward);
    R2_std = cellfun(@(R2) R2{2}, R2_trans_forward);
    hold on;
    b = bar(R2_avg, 'blue');
    errorbar(R2_avg, R2_std, 'k', 'linestyle', 'none');
    set(gca,'XTick', b.XData, 'XTickLabel', forward_names);
    title('Forward Translation');
    ylabel('R^2');
    ylim([0,1]);
    yticks(0:0.2:1);
    grid on;
    hold off;
    
    % Reverse
    subplot(2, 1, 2);
    R2_avg = cellfun(@(R2) R2{1}, R2_trans_reverse);
    R2_std = cellfun(@(R2) R2{2}, R2_trans_reverse);
    hold on;
    b = bar(R2_avg, 'blue');
    errorbar(R2_avg, R2_std, 'k', 'linestyle', 'none');
    set(gca,'XTick', b.XData, 'XTickLabel', reverse_names);
    title('Reverse Translation');
    ylabel('R^2');
    ylim([0,1]);
    yticks(0:0.2:1);
    grid on;
    hold off;
    
    saveFigures(figures, savePath);
    fprintf("%.3fs\n", toc(runTime));
        
end