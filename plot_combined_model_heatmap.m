%% Plot Saved Combined-Model Heatmap Results
clear;
close all;
clc;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Results file
results_file = 'combined_model_heatmap_results.mat';

if ~isfile(results_file)
    error('Results file not found: %s',results_file);
end

saved_results = load(results_file);

required_variables = { ...
    'K_range', ...
    'w_range', ...
    'grazer_amplitude', ...
    'cycle_period_days', ...
    'mean_stage_duration_days'};

for variable_index = 1:numel(required_variables)
    variable_name = required_variables{variable_index};

    if ~isfield(saved_results,variable_name)
        error('The MAT file does not contain %s.',variable_name);
    end
end

K_range = saved_results.K_range;
w_range = saved_results.w_range;
grazer_amplitude = saved_results.grazer_amplitude;
cycle_period_days = saved_results.cycle_period_days;
mean_stage_duration_days = ...
    saved_results.mean_stage_duration_days;

%% Check matrix dimensions
expected_size = [numel(w_range),numel(K_range)];

if ~isequal(size(grazer_amplitude),expected_size) ...
        || ~isequal(size(cycle_period_days),expected_size) ...
        || ~isequal(size(mean_stage_duration_days),expected_size)
    error(['The heatmap matrices must have size ', ...
        'length(w_range)-by-length(K_range).']);
end

%% Shared colormap
% Zero-valued non-cycles are mapped to dark blue. NaN values, including
% extinction points, are transparent and reveal the white axes background.
number_of_colors = 256;
heat_colors = parula(number_of_colors);
heat_colors(1,:) = [0.05,0.15,0.65];

%% Figure 1: grazer limit-cycle amplitude
figure('Color','w','Name','Grazer amplitude heatmap');

plot_raw_heatmap( ...
    K_range,w_range,grazer_amplitude,heat_colors, ...
    'Grazer limit-cycle amplitude', ...
    'Grazer amplitude (mg C L^{-1})');
axis square

%% Figure 2: cycle period
figure('Color','w','Name','Cycle period heatmap');

duration_edges = [0, 5, 10, 15, 20, 30, inf];

plot_binned_heatmap( ...
    K_range, ...
    w_range, ...
    cycle_period_days, ...
    duration_edges, ...
    'Mean juvenile stage duration', ...
    'Mean stage duration (days)');

%% Figure 3: mean juvenile stage duration
figure('Color','w','Name','Stage-duration heatmap');

plot_raw_heatmap( ...
    K_range,w_range,mean_stage_duration_days,heat_colors, ...
    'Mean juvenile stage duration', ...
    'Mean stage duration (days)');
axis square

%% Figure 4: simulated amplitude in (stage duration,K) space

[K_grid, ~] = meshgrid(K_range, w_range);

valid = ...
    isfinite(mean_stage_duration_days) ...
    & isfinite(grazer_amplitude);

stage_duration_values = ...
    mean_stage_duration_days(valid);

K_values = K_grid(valid);

amplitude_values = grazer_amplitude(valid);

figure( ...
    'Color', 'w', ...
    'Name', 'Stage-duration amplitude results');

scatter_handle = scatter( ...
    stage_duration_values, ...
    K_values, ...
    110, ...                 % Square-marker size
    amplitude_values, ...    % Marker color
    's', ...                 % Square markers
    'filled');

scatter_handle.MarkerEdgeColor = ...
    [0.25, 0.25, 0.25];

scatter_handle.LineWidth = 0.25;

heat_colors = parula(256);
heat_colors(1,:) = [0.05, 0.15, 0.65];

colormap(gca, heat_colors);
clim(finite_color_limits(grazer_amplitude));

amplitude_colorbar = colorbar;
amplitude_colorbar.Label.String = ...
    'Grazer amplitude (mg C L^{-1})';
amplitude_colorbar.Label.Interpreter = 'tex';

xlabel('$\overline{\tau}\;(\mathrm{days})$');
ylabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');

title('Grazer amplitude by mean stage duration');

axis square
box on

%% Local functions

function plot_raw_heatmap( ...
        K_range,w_range,data,heat_colors,plot_title,colorbar_label)

    heatmap_image = imagesc(K_range,w_range,data);

    set(heatmap_image,'AlphaData',isfinite(data));
    set(gca,'Color','w');

    axis xy
    axis square
    colormap(gca,heat_colors);
    clim(finite_color_limits(data));

    heatmap_colorbar = colorbar;
    heatmap_colorbar.Label.String = colorbar_label;
    heatmap_colorbar.Label.Interpreter = 'tex';

    xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
    ylabel('$w$');
    title(plot_title);

    box on
end

function plot_binned_heatmap( ...
        x_values, y_values, data, bin_edges, ...
        plot_title, colorbar_label)

    % Convert continuous values into numbered bins
    binned_data = discretize(data, bin_edges);

    number_of_bins = numel(bin_edges)-1;

    heatmap_image = imagesc( ...
        x_values, ...
        y_values, ...
        binned_data);

    % Missing/extinction values remain white
    set(heatmap_image, ...
        'AlphaData', isfinite(binned_data));

    set(gca, 'Color', 'w');

    axis xy
    axis square

    % One distinct color per bin
    bin_colors = parula(number_of_bins);
    colormap(gca, bin_colors);

    clim([0.5, number_of_bins+0.5]);

    %% Generate colorbar labels automatically

    bin_labels = cell(1, number_of_bins);

    for bin_index = 1:number_of_bins

        lower_edge = bin_edges(bin_index);
        upper_edge = bin_edges(bin_index+1);

        if isinf(upper_edge)
            bin_labels{bin_index} = sprintf( ...
                '$>%g$', lower_edge);
        else
            bin_labels{bin_index} = sprintf( ...
                '$%g$--$%g$', ...
                lower_edge, upper_edge);
        end
    end

    heatmap_colorbar = colorbar;
    heatmap_colorbar.Ticks = 1:number_of_bins;
    heatmap_colorbar.TickLabels = bin_labels;
    heatmap_colorbar.TickLabelInterpreter = 'latex';

    heatmap_colorbar.Label.String = colorbar_label;
    heatmap_colorbar.Label.Interpreter = 'tex';

    xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
    ylabel('$w$');
    title(plot_title);

    box on
end

function limits = finite_color_limits(data)
    finite_data = data(isfinite(data));

    if isempty(finite_data)
        limits = [0,1];
        return
    end

    lower_limit = min(0,min(finite_data));
    upper_limit = max(finite_data);

    if upper_limit <= lower_limit
        upper_limit = lower_limit+1;
    end

    limits = [lower_limit,upper_limit];
end
