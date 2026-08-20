%% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_matthiessen_ZrNiSn_kScan_holes.mat';

%% Initialize for plotting
figure;
fig = gcf;
fig.Position(3:4) = [800, 600];
hold on;

numBands = 3;  % Adjust this to your actual number of bands
palette = [
    0.1216 0.4667 0.7059;  % blue
    1.0000 0.4980 0.0549;  % orange
    0.1725 0.6275 0.1725;  % green
    0.8392 0.1529 0.1569;  % red
    0.5804 0.4039 0.7412;  % purple
    0.5490 0.3373 0.2941;  % brown
    0.8902 0.4667 0.7608;  % pink
    0.4980 0.4980 0.4980;  % gray
];

h = gobjects(numBands, 1);
h_created = false(numBands, 1);  % Track which handles were actually created

% Initialize arrays for data collection
E_avg_all = [];          % For average line
sca_avg_all = [];
all_E = [];              % For individual k-points
all_sca = [];            % For individual k-points
all_band = [];           % Track band index for each point

% Configuration
targetTempIdx = 1;       % Temperature index (1=300K, 2=400K, etc.)
tempVal = 300 + (targetTempIdx-1)*100;

fprintf('Processing Matthiessen data for T = %d K\n', tempVal);

%% Main processing loop
% Find the maximum index based on your data structure
% For Matthiessen data, it's usually organized by state_ID length
max_idx = length(state_ID);

for i = 1:max_idx
    E_row = [];
    sca_row = [];
    
    for b = 1:numBands
        % Check if this band exists for this index
        if size(state_ID, 2) >= b && i <= size(state_ID, 1)
            % Check if data exists for this state and band
            if ~isempty(taus_matth(i).x)  % Adjust this based on your structure
                
                % Get energy for this point
                E = state_ID(i).E;
                
                % Extract tau for the specific temperature
                % For Matthiessen, the structure might be different
                % Try these options:
                
                % Option 1: If taus_matth is similar to taus_POP
                if ndims(taus_matth(i).x) == 3
                    % 3D array: (something, temperature, k-points)
                    tau_k = squeeze(taus_matth(i).x(8, targetTempIdx, :));
                elseif ndims(taus_matth(i).x) == 2
                    % 2D array: (k-points, temperature)
                    if size(taus_matth(i).x, 2) >= targetTempIdx
                        tau_k = taus_matth(i).x(:, targetTempIdx);
                    else
                        tau_k = taus_matth(i).x(:, 1);  % Default to first column
                    end
                else
                    % 1D array: just k-points
                    tau_k = taus_matth(i).x(:);
                end
                
                % Filter: Keep only positive tau values
                positive_mask = tau_k > 0;
                
                if any(positive_mask)
                    tau_k_positive = tau_k(positive_mask);
                    sca_k = 1 ./ tau_k_positive;  % Scattering rates
                    
                    % Store individual k-points for saving
                    n_points = length(sca_k);
                    all_E = [all_E; repmat(E, n_points, 1)];
                    all_sca = [all_sca; sca_k];
                    all_band = [all_band; repmat(b, n_points, 1)];
                    
                    % Scatter plot for all k-points at this energy
                    htmp = scatter(repmat(E, size(sca_k)), sca_k, 30, palette(b,:), 'filled');
                    
                    % Average over k-points for this band at this energy
                    sca_band_avg = mean(sca_k);
                    scatter(E, sca_band_avg, 50, palette(b,:), 'filled', 'd');
                    
                    % Store handle only if this is the first valid data for this band
                    if ~h_created(b)
                        h(b) = htmp;
                        h_created(b) = true;
                    end
                    
                    % Collect for total average over bands
                    E_row(end+1) = E;
                    sca_row(end+1) = sca_band_avg;
                end
            end
        end
    end
    
    % After looping over bands, take the overall average
    if ~isempty(E_row)
        E_avg_all(end+1) = mean(E_row);
        sca_avg_all(end+1) = mean(sca_row);
    end
    
    % Progress indicator
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', i, max_idx);
    end
end

%% Overlay average line
avg_line = plot(E_avg_all, sca_avg_all, 'k-', 'LineWidth', 2.5);

%% Format the plot
set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([0, 0.45]);
grid on;

xlabel('Energy [eV]', 'FontSize', 12);
ylabel('\tau^{-1} [s^{-1}]', 'FontSize', 12);
title(sprintf('ZrNiSn Matthiessen - T=%dK, all bands', tempVal), 'FontSize', 14);

%% Create legend - FIXED VERSION
legend_handles = [];
legend_labels = {};

% Add band handles that were actually created
for b = 1:numBands
    if h_created(b)
        legend_handles = [legend_handles; h(b)];
        legend_labels{end+1} = sprintf('Band %d', b);
    end
end

% Add average line
legend_handles = [legend_handles; avg_line];
legend_labels{end+1} = 'Average';

% Only create legend if we have valid handles
if ~isempty(legend_handles)
    lgd = legend(legend_handles, legend_labels);
    lgd.FontSize = 10;
    lgd.Location = 'best';
else
    warning('No valid data to create legend');
end

%% Save data for xmgrace plotting
fprintf('\n=== Saving data for xmgrace ===\n');
baseName = sprintf('ZrNiSn_Matthiessen_%dK', tempVal);

% Create output directory
outputDir = sprintf('matthiessen_xmgrace_%dK', tempVal);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 1. Save all individual k-points (for scatter plots)
if ~isempty(all_E)
    individual_data = [all_E, all_sca];
    filename = fullfile(outputDir, sprintf('%s_individual.dat', baseName));
    dlmwrite(filename, individual_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Saved %d individual k-points to %s\n', length(all_E), filename);
end

% 2. Save energy-averaged data (for line plots)
if ~isempty(E_avg_all)
    avg_data = [E_avg_all(:), sca_avg_all(:)];
    filename = fullfile(outputDir, sprintf('%s_average.dat', baseName));
    dlmwrite(filename, avg_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Saved %d averaged points to %s\n', length(E_avg_all), filename);
end

% 3. Save with band information (3 columns: E, sca, band)
if ~isempty(all_E)
    band_data = [all_E, all_sca, all_band];
    filename = fullfile(outputDir, sprintf('%s_with_bands.dat', baseName));
    dlmwrite(filename, band_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Saved band-labeled data to %s\n', filename);
end

% 4. Save formatted file with headers
if ~isempty(E_avg_all)
    filename = fullfile(outputDir, sprintf('%s_formatted.dat', baseName));
    fid = fopen(filename, 'w');
    fprintf(fid, '# ZrNiSn Matthiessen Scattering Rates\n');
    fprintf(fid, '# Temperature: %d K\n', tempVal);
    fprintf(fid, '# Energy (eV)  Scattering Rate (s^-1)\n');
    for i = 1:length(E_avg_all)
        fprintf(fid, '%.6f  %.6e\n', E_avg_all(i), sca_avg_all(i));
    end
    fclose(fid);
    fprintf('✓ Saved formatted data to %s\n', filename);
end

% 5. Save per-band average files
fprintf('\n--- Per-band average files ---\n');
for b = 1:numBands
    if h_created(b)
        % Extract data for this specific band
        band_mask = (all_band == b);
        E_band = all_E(band_mask);
        sca_band = all_sca(band_mask);
        
        if ~isempty(E_band)
            % Get unique energies for this band
            [unique_E_band, ~, idx_E] = unique(E_band);
            E_band_avg = unique_E_band;
            sca_band_avg = accumarray(idx_E, sca_band, [], @mean);
            
            filename = fullfile(outputDir, sprintf('%s_band%d_average.dat', baseName, b));
            fid = fopen(filename, 'w');
            fprintf(fid, '# Band %d Average Scattering Rates\n', b);
            fprintf(fid, '# Energy (eV)  Scattering Rate (s^-1)\n');
            for j = 1:length(E_band_avg)
                fprintf(fid, '%.6f  %.6e\n', E_band_avg(j), sca_band_avg(j));
            end
            fclose(fid);
            fprintf('✓ Band %d: %d points saved to %s\n', b, length(E_band_avg), filename);
        end
    end
end

% 6. Create xmgrace batch file
filename = fullfile(outputDir, sprintf('%s_plot.bfile', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# xmgrace batch file for Matthiessen scattering rates\n');
fprintf(fid, 'title "ZrNiSn Matthiessen - T = %dK"\n', tempVal);
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [s\\S-1\\N]"\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'legend on\n\n');

% Add per-band averages
set_count = 0;
for b = 1:numBands
    if h_created(b)
        fprintf(fid, 'READ NXY "band%d_average.dat"\n', b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d legend "Band %d"\n', set_count, b);
        set_count = set_count + 1;
    end
end

% Add overall average
fprintf(fid, 'READ NXY "average.dat"\n');
fprintf(fid, 's%d line color 1\n', set_count);
fprintf(fid, 's%d line linewidth 3\n', set_count);
fprintf(fid, 's%d legend "Overall Average"\n', set_count);

fprintf(fid, '\nAUTOSCALE\n');
fclose(fid);
fprintf('✓ xmgrace batch file saved to %s\n', filename);

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('Temperature: %d K\n', tempVal);
fprintf('Processed %d energy points\n', max_idx);
fprintf('Found valid data in %d bands\n', sum(h_created));
if ~isempty(all_E)
    fprintf('Total individual k-points: %d\n', length(all_E));
end
if ~isempty(E_avg_all)
    fprintf('Average curve points: %d\n', length(E_avg_all));
end
fprintf('All files saved in: %s\n', outputDir);
fprintf('To plot in xmgrace: cd %s && xmgrace -batch %s_plot.bfile\n', outputDir, baseName);