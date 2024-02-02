function plot_pattern_slice(yvector, xvector, ...
                           yaxislabel, xaxislabel, ...
                           yaxismin, yaxismax, ...
                           xaxismin, xaxismax, color)
    % plot_pattern_slice - Plots a 2D pattern slice
    %
    %   plot_pattern_slice(yvector, xvector, yaxislabel, xaxislabel,
    %                      yaxismin, yaxismax, xaxismin, xaxismax, color)
    %
    %   Plots a 2D pattern slice represented by the vectors yvector and
    %   xvector, with specified axis labels, axis limits, and color.
    %
    % Input:
    %   yvector - Vector representing the y-values of the pattern
    %   xvector - Vector representing the x-values of the pattern
    %   yaxislabel - Label for the y-axis
    %   xaxislabel - Label for the x-axis
    %   yaxismin - Minimum value for the y-axis
    %   yaxismax - Maximum value for the y-axis
    %   xaxismin - Minimum value for the x-axis
    %   xaxismax - Maximum value for the x-axis
    %   color - Color specification for the plot
    %
    % Output:
    %   The function produces a 2D plot of the pattern slice.
    %

    % Plot the pattern slice
    plot(xvector, yvector, color, 'LineWidth', 2)

    % Add grid to the plot
    grid on

    % Set axis labels
    xlabel(xaxislabel)
    ylabel(yaxislabel)

    % Set axis limits
    xlim([xaxismin, xaxismax])
    ylim([yaxismin, yaxismax])
end
