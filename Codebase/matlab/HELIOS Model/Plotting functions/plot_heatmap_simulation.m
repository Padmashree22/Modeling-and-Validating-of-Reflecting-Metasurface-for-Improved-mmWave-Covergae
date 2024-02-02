function plot_heatmap_simulation(matrix, xvector, yvector, ...
                                 colorbarlabel, xaxislabel, yaxislabel, ...
                                 colorbarmin, colorbarmax)
    % plot_heatmap_simulation - Plot a heatmap for simulation results
    %
    %   plot_heatmap_simulation(matrix, xvector, yvector, colorbarlabel,
    %                            xaxislabel, yaxislabel, colorbarmin,
    %                            colorbarmax) generates a heatmap plot for
    %   simulation results given a matrix of values, x and y vectors,
    %   colorbar label, x-axis label, y-axis label, colorbar limits.
    %
    % Input:
    %   matrix - RCS Matrix
    %   xvector - Vector representing reflected wave
    %   yvector - Vector representing reflected wave
    %   colorbarlabel - Label for the colorbar
    %   xaxislabel - Label for the x-axis
    %   yaxislabel - Label for the y-axis
    %   colorbarmin - Minimum value for the colorbar
    %   colorbarmax - Maximum value for the colorbar
    %

    % Create heatmap using imagesc, flipud, and set YDir
    im = imagesc(flipud(matrix'), [colorbarmin, colorbarmax]);
    set(gca, 'YDir', 'normal')

    % Add colorbar and set colormap to jet
    c = colorbar;
    colormap jet;

    % Set colormap explicitly for the figure
    map = colormap;
    set(gcf, 'Colormap', map);

    % Label the colorbar
    ylabel(c, colorbarlabel);

    % Set X and Y limits for the heatmap
    im.Parent.XLim = [min(xvector), max(xvector)];
    im.Parent.YLim = [min(yvector), max(yvector)];
    im.XData = [min(xvector), max(xvector)];
    im.YData = [min(yvector), max(yvector)];

    % Label the axes
    xlabel(xaxislabel);
    ylabel(yaxislabel);

    % Turn off legend box
    legend('boxoff');
end