function plot_heatmap_analytical(matrix, xvector, yvector, ...
                                 colorbarlabel, xaxislabel, yaxislabel, ...
                                 colorbarmin, colorbarmax)
    % plot_heatmap_analytical - Create a heatmap plot for analytical model
    %
    %   plot_heatmap_analytical(matrix, xvector, yvector, ...
    %                            colorbarlabel, xaxislabel, yaxislabel, ...
    %                            colorbarmin, colorbarmax)
    %
    %   This function creates a heatmap plot from a matrix with specified
    %   analytical options, including color mapping and handling of NaN
    %   values.
    %
    %   Parameters:
     %   matrix - RCS Matrix
    %   xvector - Vector representing reflected wave
    %   yvector - Vector representing reflected wave
    %   - colorbarlabel: Label for the colorbar
    %   - xaxislabel: Label for the x-axis
    %   - yaxislabel: Label for the y-axis
    %   - colorbarmin: Minimum value for the colorbar
    %   - colorbarmax: Maximum value for the colorbar
    %

    nanColor = [0, 0, 0];           % Black color for NaN
    greyColor = [0.5, 0.5, 0.5];    % Grey color for the other region
    map = colormap(jet(256));       % Use the colormap you prefer
  
    nanColorIndex = size(map, 1) + 1;  % Index for NaN color
    infColorIndex = nanColorIndex + 1;  % Index for grey color
   
    map = [nanColor;  map; greyColor;];

    % Create the heatmap
    im = imagesc(matrix, [colorbarmin, colorbarmax]);

    % Set NaN color and handle inf values
    cdata = im.CData;
    cdata(isnan(cdata)) = nanColorIndex;
    cdata(cdata < colorbarmin & ~isinf(cdata)  & isinf(cdata)) = infColorIndex;
    im.CData = cdata;
    
    % Set axis direction, colormap, and colorbar
    set(gca, 'YDir', 'normal')
    c = colorbar;
    set(gcf, 'Colormap', map);
    
    % Set colorbar label and axis limits
    ylabel(c, colorbarlabel);
    im.Parent.XLim = [min(xvector), max(xvector)];
    im.Parent.YLim = [min(yvector), max(yvector)];
    im.XData = [min(xvector), max(xvector)];
    im.YData = [min(yvector), max(yvector)];
    
    % Set axis labels and turn off legend box
    xlabel(xaxislabel);
    ylabel(yaxislabel);
    legend('boxoff');
end
