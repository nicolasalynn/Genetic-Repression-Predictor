function data_plotter(X, y_actual, y_predicted, m)


    together = [X; y_actual; y_predicted]';
    together = sortrows(together);


    real_mean = mean(y_actual);
    pred_mean = mean(y_predicted);
    range_vals = range(X);
    correlation_y = (corr(together(:, 2), together(:, 3), 'type', 'pearson'))*100;

    fprintf("\n######################################################\n")
    fprintf("Slope of Regression:\t\t%.2d\nMean of Actual Data: \t\t%.2d\nMean of Predicted Data: " + ...
        "\t%d\nRange of X Values: \t\t%d\n" +...
        "Correlation of Y Values:\t%2.2f%%", ...
        m, real_mean, pred_mean, range_vals(1), correlation_y);
    fprintf("\n######################################################\n\n")

    unique_vals = unique(together(:, 1));
    
    if length(unique_vals) > 10
        
        figure(1)
        plot(X, y_predicted, 'r')
        hold on
        scatter(X(1:250), y_actual(1:250), 'b')
        hold off

    end

    if length(unique_vals) < 3
        unique_table = zeros(length(unique_vals), 4);
        for i = 1:length(unique_vals)
            unique_table(i, 1) = unique_vals(i);
            for j = 1:size(together, 1)
                if together(j, 1) == unique_vals(i)
                    unique_table(i, 2) = unique_table(i, 2) + together(j, 2);
                    unique_table(i, 3) = unique_table(i, 3) + together(j, 3);
                    unique_table(i, 4) = unique_table(i, 4) + 1;
                end
            end
        end

        average_actual_per_val = unique_table(:, 2) ./ unique_table(:, 4);
        average_pred_per_val = unique_table(:, 3) ./ unique_table(:, 4);
        average_per_val = [unique_table(:, 1), average_actual_per_val, average_pred_per_val];

        fprintf("#####################################################\n");
        fprintf("X Value \tMean Predicted \t\tMean Actual\n");
        for k = 1:size(average_per_val, 1)
            fprintf("%i\t\t%d\t\t\t%d\n", average_per_val(k, 1), average_per_val(k, 2), average_per_val(k, 3))
        end
        fprintf("#####################################################\n");

        plot(average_per_val(:, 1), average_per_val(:, 2))
        hold on
        plot(average_per_val(:, 1), average_per_val(:, 3))
        hold off
    end

    
end