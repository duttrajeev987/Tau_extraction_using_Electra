%% 1. Load the data files
load 'TE_ZrNiSn_kScan_holes.mat'; 
load 'RelaxTimes_matthiessen_ZrNiSn_kScan_holes.mat'; 

%% 2. Configuration & Parameters
targetTempIdx = 1;      % 1=300K, 2=400K, 3=500K, 4=600K, etc.
threshold = 1e15;       % Threshold: Ignore rates > 10^15 for the average
tempVal = 300 + (targetTempIdx-1)*100;

% Manual band configuration
numBands = 3;  % Set this to your actual number of bands
% Define band indices (if bands are numbered differently, adjust here)
band_indices = 1:numBands;  % or manually: [1, 2, 3, 4, 5, 6, 7, 8]

% Color palette for multiple bands
palette = [
    0.1216 0.4667 0.7059;  % blue
    1.0000 0.4980 0.0549;  % orange
    0.1725 0.6275 0.1725;  % green
    0.8392 0.1529 0.1569;  % red
    0.5804 0.4039 0.7412;  % purple
    0.5490 0.3373 0.2941;  % brown
    0.8902 0.4667 0.7608;  % pink
    0.4980 0.4980 0.4980;  % gray
    0.7373 0.7412 0.1333;  % olive (extra if needed)
    0.0902 0.7451 0.8118;  % cyan (extra if needed)
];

% Visualization Setup
figure;
fig = gcf;
fig.Position(3:4) = [800, 600];
hold on;

% Initialize arrays for legend (using simpler approach)
plot_handles = [];
legend_strings = {};

% Initialize data containers for ALL bands combined
all_E_combined   = [];          
all_sca_combined = [];

% Initialize cell arrays to store per-band data
E_avg_by_band = cell(numBands, 1);
sca_avg_by_band = cell(numBands, 1);
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);

%% 3. Main Processing Loop with Manual Band Separation
fprintf('\n=== Processing %d bands at %dK ===\n', numBands, tempVal);

for bandIdx = 1:numBands
    current_band = band_indices(bandIdx);
    fprintf('Processing band %d/%d (band index %d)...\n', bandIdx, numBands, current_band);
    
    % Initialize containers for this band
    E_avg_band = [];
    sca_avg_band = [];
    E_all_band = [];
    sca_all_band = [];
    
    % Iterating through all energy slices
    for i = 1:length(taus_matth)
        
        % Check if this slice belongs to current band
        % Adjust this logic based on your actual data structure
        if isfield(state_ID, 'band')
            % If state_ID has a 'band' field
            if state_ID(i).band ~= current_band
                continue;
            end
        else
            % Alternative: if bands are encoded in state_ID differently
            % For example, if stored in a separate array or field
            % You might need to modify this section
            continue;  % Skip if no band information available
        end
        
        % Check if data exists for this energy slice and temperature column
        if ~isempty(taus_matth(i).x) && size(taus_matth(i).x, 2) >= targetTempIdx
            
            % Get energy for this slice
            E = state_ID(i).E;
            
            % Extract tau for the specific temperature (all k-points in rows)
            tau_k = taus_matth(i).x(:, targetTempIdx);
            
            % Filter 1: Keep only positive tau values
            positive_mask = tau_k > 0;
            
            if any(positive_mask)
                tau_k_positive = tau_k(positive_mask);
                sca_k = 1 ./ tau_k_positive;  % Scattering rates (1/tau)
                
                % Filter 2: Apply threshold for average calculation
                avg_mask = sca_k <= threshold;
                sca_k_for_avg = sca_k(avg_mask);
                
                % Collect ALL individual points for this band
                E_all_band = [E_all_band; repmat(E, length(sca_k), 1)];
                sca_all_band = [sca_all_band; sca_k];
                
                % Add to combined dataset
                all_E_combined = [all_E_combined; repmat(E, length(sca_k), 1)];
                all_sca_combined = [all_sca_combined; sca_k];
                
                % Calculate average for this band
                if ~isempty(sca_k_for_avg)
                    sca_slice_avg = mean(sca_k_for_avg);
                    
                    E_avg_band(end+1) = E;
                    sca_avg_band(end+1) = sca_slice_avg;
                end
            end
        end
    end
    
    % Store band-specific data
    E_avg_by_band{bandIdx} = E_avg_band;
    sca_avg_by_band{bandIdx} = sca_avg_band;
    all_E_by_band{bandIdx} = E_all_band;
    all_sca_by_band{bandIdx} = sca_all_band;
    
    % Plot band average (bold line) and collect handle for legend
    if ~isempty(E_avg_band)
        % Plot individual k-points as light scatter (optional)
        % scatter(E_all_band, sca_all_band, 15, palette(bandIdx,:), 'filled', 'MarkerFaceAlpha', 0.05);
        
        % Plot the average line and store the handle
        h = plot(E_avg_band, sca_avg_band, '-', 'Color', palette(bandIdx,:), ...
                'LineWidth', 2.5, 'DisplayName', sprintf('Band %d', current_band));
        
        % Store handle and legend string
        plot_handles(end+1) = h;
        legend_strings{end+1} = sprintf('Band %d', current_band);
        
        fprintf('  Band %d: %d energy points averaged\n', current_band, length(E_avg_band));
    end
end

%% 4. Calculate and Plot Overall Average (across all bands)
fprintf('\nCalculating overall average across all bands...\n');

% Group all data by energy and calculate overall average
[unique_E, ~, E_groups] = unique(all_E_combined);
overall_E_avg = [];
overall_sca_avg = [];

for i = 1:length(unique_E)
    % Get all scattering rates at this energy across all bands
    sca_at_E = all_sca_combined(E_groups == i);
    
    % Apply threshold for overall average
    sca_at_E_filtered = sca_at_E(sca_at_E <= threshold);
    
    if ~isempty(sca_at_E_filtered)
        overall_E_avg(end+1) = unique_E(i);
        overall_sca_avg(end+1) = mean(sca_at_E_filtered);
    end
end

% Plot overall average (thick black line)
if ~isempty(overall_E_avg)
    h_overall = plot(overall_E_avg, overall_sca_avg, 'k-', 'LineWidth', 3.5, ...
                     'DisplayName', 'Overall Average');
    plot_handles(end+1) = h_overall;
    legend_strings{end+1} = 'Overall Average';
end

%% 5. Finalize Plot
set(gca, 'YScale', 'log');
grid on;
ylim([1e11, 2e15]);
xlim([0, 0.45]);
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Scattering Rate \tau^{-1} [s^{-1}]', 'FontSize', 12);
title(sprintf('ZrNiSn Matthiessen - %d Bands | T = %dK | Filter: < 10^{15}', ...
      numBands, tempVal), 'FontSize', 14);

% Create legend - FIXED VERSION
if ~isempty(plot_handles)
    % Make sure all handles are valid graphics objects
    valid_idx = isgraphics(plot_handles);
    if any(valid_idx)
        legend(plot_handles(valid_idx), legend_strings(valid_idx), ...
               'Location', 'best', 'FontSize', 10);
    else
        warning('No valid plot handles for legend');
    end
end

%% 6. Export Data for xmgrace - Organized by Bands
fprintf('\n=== Exporting Data for xmgrace ===\n');
baseName = sprintf('ZrNiSn_Matth_%dK', tempVal);

% Create output directory
outputDir = sprintf('xmgrace_data_%dK', tempVal);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% A. Save Individual band data files
fprintf('\n--- Band-specific files ---\n');
for bandIdx = 1:numBands
    current_band = band_indices(bandIdx);
    
    % Individual k-points for this band
    if ~isempty(all_E_by_band{bandIdx})
        filename_ind = fullfile(outputDir, sprintf('%s_band%02d_individual.dat', baseName, current_band));
        dlmwrite(filename_ind, [all_E_by_band{bandIdx}, all_sca_by_band{bandIdx}], ...
                'delimiter', '\t', 'precision', '%.6e');
        fprintf('  ✓ Band %d individual: %s (%d points)\n', current_band, filename_ind, ...
                length(all_E_by_band{bandIdx}));
    end
    
    % Filtered averages for this band (the trend line)
    if ~isempty(E_avg_by_band{bandIdx})
        filename_avg = fullfile(outputDir, sprintf('%s_band%02d_average.dat', baseName, current_band));
        
        % Create formatted file with header
        fid = fopen(filename_avg, 'w');
        fprintf(fid, '# ZrNiSn Band %d Scattering Rate Average\n', current_band);
        fprintf(fid, '# Temperature: %d K\n', tempVal);
        fprintf(fid, '# Filter: rates > %.1e excluded\n', threshold);
        fprintf(fid, '# Energy(eV)  ScatteringRate(s^-1)\n');
        for i = 1:length(E_avg_by_band{bandIdx})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band{bandIdx}(i), sca_avg_by_band{bandIdx}(i));
        end
        fclose(fid);
        
        fprintf('  ✓ Band %d average: %s (%d points)\n', current_band, filename_avg, ...
                length(E_avg_by_band{bandIdx}));
    end
end

% B. Save Overall Average (across all bands)
fprintf('\n--- Combined files ---\n');
if ~isempty(overall_E_avg)
    filename_overall = fullfile(outputDir, sprintf('%s_overall_average.dat', baseName));
    
    fid = fopen(filename_overall, 'w');
    fprintf(fid, '# ZrNiSn Overall Scattering Rate Average (All Bands)\n');
    fprintf(fid, '# Temperature: %d K\n', tempVal);
    fprintf(fid, '# Filter: rates > %.1e excluded\n', threshold);
    fprintf(fid, '# Energy(eV)  ScatteringRate(s^-1)\n');
    for i = 1:length(overall_E_avg)
        fprintf(fid, '%.6f  %.6e\n', overall_E_avg(i), overall_sca_avg(i));
    end
    fclose(fid);
    
    fprintf('  ✓ Overall average: %s (%d points)\n', filename_overall, length(overall_E_avg));
end

% C. Create xmgrace batch file for easy plotting
filename_batch = fullfile(outputDir, sprintf('%s_multi_band.bfile', baseName));
fid = fopen(filename_batch, 'w');

fprintf(fid, '# xmgrace batch file for multi-band plot\n');
fprintf(fid, '# Generated automatically from MATLAB\n\n');
fprintf(fid, 'WORLD XMIN 0\n');
fprintf(fid, 'WORLD XMAX 0.45\n');
fprintf(fid, 'WORLD YMIN 1e11\n');
fprintf(fid, 'WORLD YMAX 2e15\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [s\\S-1\\N]"\n');
fprintf(fid, 'title "ZrNiSn Matthiessen - %d Bands, T = %dK"\n', numBands, tempVal);
fprintf(fid, 'subtitle "Filter: rates < %.1e s\\S-1\\N"\n', threshold);
fprintf(fid, 'legend on\n\n');

% Add band data
set_count = 0;
for bandIdx = 1:numBands
    current_band = band_indices(bandIdx);
    if ~isempty(E_avg_by_band{bandIdx})
        fprintf(fid, 'READ NXY "band%02d_average.dat"\n', current_band);
        fprintf(fid, 's%d line color %d\n', set_count, bandIdx);
        fprintf(fid, 's%d line linewidth 2.5\n', set_count);
        fprintf(fid, 's%d legend "Band %d"\n', set_count, current_band);
        set_count = set_count + 1;
    end
end

% Add overall average
if ~isempty(overall_E_avg)
    fprintf(fid, 'READ NXY "overall_average.dat"\n');
    fprintf(fid, 's%d line color 1\n', set_count);  % Black
    fprintf(fid, 's%d line linewidth 3.5\n', set_count);
    fprintf(fid, 's%d legend "Overall Average"\n', set_count);
end

fprintf(fid, '\nAUTOSCALE\n');
fclose(fid);
fprintf('  ✓ xmgrace batch file: %s\n', filename_batch);

% D. Create summary file
filename_summary = fullfile(outputDir, sprintf('%s_README.txt', baseName));
fid = fopen(filename_summary, 'w');
fprintf(fid, '========================================\n');
fprintf(fid, 'ZrNiSn Scattering Rate Data Summary\n');
fprintf(fid, '========================================\n');
fprintf(fid, 'Temperature: %d K\n', tempVal);
fprintf(fid, 'Filter threshold: %.1e s^-1\n', threshold);
fprintf(fid, 'Number of bands: %d\n', numBands);
fprintf(fid, 'Band indices: %s\n', mat2str(band_indices));
fprintf(fid, '\n--- File Description ---\n');
fprintf(fid, '*_bandXX_individual.dat : All k-point scattering rates for band XX\n');
fprintf(fid, '*_bandXX_average.dat    : Energy-averaged scattering rates for band XX\n');
fprintf(fid, '*_overall_average.dat   : Average across all bands\n');
fprintf(fid, '*_multi_band.bfile      : xmgrace batch file for plotting\n');
fprintf(fid, '\n--- Band Statistics ---\n');
for bandIdx = 1:numBands
    current_band = band_indices(bandIdx);
    if ~isempty(all_sca_by_band{bandIdx})
        fprintf(fid, 'Band %2d: %6d individual points, %4d energy averages\n', ...
                current_band, length(all_sca_by_band{bandIdx}), ...
                length(sca_avg_by_band{bandIdx}));
    else
        fprintf(fid, 'Band %2d: No data\n', current_band);
    end
end
fprintf(fid, 'Total combined points: %d\n', length(all_sca_combined));
fclose(fid);
fprintf('  ✓ Summary file: %s\n', filename_summary);

%% 7. Final Summary
fprintf('\n========================================\n');
fprintf('Processing Complete!\n');
fprintf('========================================\n');
fprintf('Output directory: %s\n', outputDir);
fprintf('Files created:\n');
fprintf('  - %d band-specific average files\n', numBands);
fprintf('  - %d band-specific individual files\n', numBands);
fprintf('  - 1 overall average file\n');
fprintf('  - 1 xmgrace batch file\n');
fprintf('  - 1 README file\n');
fprintf('\nTo plot in xmgrace, use:\n');
fprintf('  xmgrace -batch %s\n', filename_batch);
fprintf('========================================\n');