
% Parameters
a = 0.5;
b_values = linspace(0, 2, 500); % Range of b values
z_max = []; % To store the maxima of z
colors = []; % To store colors for each point

% Initial conditions
x0 = [0.4; 0.5; 0.5]; % Initial condition

% Create colormap
cmap = jet(256); % Using jet colormap for vibrant colors

% Loop over the range of b values
for i = 1:length(b_values)
    b = b_values(i);
    
    % Define the system of ODEs
    system = @(t, x) [x(2) - x(1);
                      -x(1) * x(3) + a * x(2);
                      -b + x(1) * x(2)];
    % Time span
    tspan = [0 300];
    
    % Solve the system
    [~, sol] = ode45(system, tspan, x0);
    
    % Consider only the second half to remove transients
    z_data = sol(round(length(sol)/2):end, 3);
    
    % Find the local maxima of z
    [peaks, locs] = findpeaks(z_data);
    num_peaks = length(peaks);
    
    % Store the results
    z_max = [z_max; b*ones(num_peaks, 1), peaks];
    
    % Assign colors based on the z-value
    min_z = 0;
    max_z = 2; % Adjust this value based on your expected z range
    normalized_z = (peaks - min_z) / (max_z - min_z);
    color_indices = max(1, min(256, round(normalized_z * 255) + 1));
    colors = [colors; cmap(color_indices, :)];
end

% Create figure with white background
figure('Color', 'white');

% Plot bifurcation diagram using scatter with colors
scatter(z_max(:, 1), z_max(:, 2), 2, colors, 'filled');

% Customize the plot
title('Bifurcation diagram of the Pehlivan system', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('b', 'FontSize', 12);
ylabel('z', 'FontSize', 12);
box on;

% % Add colorbar with correct tick labels
% colorbar('Ticks', linspace(0, 1, 5), ...
%          'TickLabels', linspace(min_z, max_z, 5));

% Make the plot more visually appealing
set(gca, 'LineWidth', 1.2, 'Box', 'on');

% ================= AUTO-SAVE SETTINGS =================
outputDir = 'Figures_Bifurcation';
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

dpi = 600;
filename = fullfile(outputDir,'Bifurcation_Pehlivan_Integer.png');

print(gcf, filename, '-dpng', ['-r',num2str(dpi)]);

fprintf('Bifurcation diagram saved as:\n%s\n', filename);


