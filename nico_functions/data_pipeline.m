function [X, y_obs, y_pred, m, correlation] = data_pipeline(X_ugly, y_ugly)

    size_x = size(X_ugly);
    size_y = size(y_ugly);
    
    if (size_x(1) == 1 || size_x(2) == 1) && (size_y(1) == 1 || size_y(2) == 1)
        X = X_ugly;
        y_obs = y_ugly;
    else
        [X, y_obs] = create_usable_data(X_ugly, y_ugly);
    end
    
    m = regress(y_obs', X');
    y_pred =  X .* m;
    
    correlation = corr(y_pred', y_obs', 'type', 'spearman');
    data_plotter(X, y_obs, y_pred, m);
    
end
