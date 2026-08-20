%% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';

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
E_avg_all = [];          % For average line (across all bands)
sca_avg_all = [];

% Initialize per-band storage
E_avg_by_band = cell(numBands, 1);    % Average per band
sca_avg_by_band = cell(numBands, 1);
all_E_by_band = cell(numBands, 1);    % Individual points per band
all_sca_by_band = cell(numBands, 1);

% For combined storage
all_E = [];              % For individual k-points (all bands)
all_sca = [];            % For individual k-points (all bands)
all_band = [];           % Track band index for each point

% Configuration
EF_index = 8;            % Index for EF=0
T_index = 7;             % Temperature index (1=900K, adjust if needed)

fprintf('Processing IIS data for EF=0, T=900K\n');

%% Main processing loop
max_idx = 451;  % Fixed 451 energy points

for i = 1:max_idx
    E_row = [];
    sca_row = [];
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            % Get energy for this point
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= EF_index && size(tau_3D, 2) >= T_index
                % Extract tau(EF=0, T=900K, all k-points)
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                
                % Filter: Keep only positive tau values
                positive_mask = tau_k > 0;
                
                if any(positive_mask)
                    tau_k_positive = tau_k(positive_mask);
                    sca_k = 1 ./ tau_k_positive;  % Scattering rates
                    
                    % Store individual k-points for this band
                    n_points = length(sca_k);
                    all_E_by_band{b} = [all_E_by_band{b}; repmat(E, n_points, 1)];
                    all_sca_by_band{b} = [all_sca_by_band{b}; sca_k];
                    
                    % Also store in combined arrays
                    all_E = [all_E; repmat(E, n_points, 1)];
                    all_sca = [all_sca; sca_k];
                    all_band = [all_band; repmat(b, n_points, 1)];
                    
                    % Scatter plot for all k-points at this energy
                    htmp = scatter(repmat(E, size(sca_k)), sca_k, 30, palette(b,:), 'filled');
                    
                    % Average over k-points for this band at this energy
                    sca_band_avg = mean(sca_k);
                    scatter(E, sca_band_avg, 50, palette(b,:), 'filled', 'd');
                    
                    % Store average for this band
                    E_avg_by_band{b}(end+1) = E;
                    sca_avg_by_band{b}(end+1) = sca_band_avg;
                    
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
ylabel('\tau^{-1} [fs^{-1}]', 'FontSize', 12);
title('ZrNiSn IIS - EF=0, T=900K, all bands', 'FontSize', 14);

%% Create legend
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

%% ===== SAVE DATA FOR XMGRACE - SEPARATE FILES PER BAND =====
fprintf('\n=== Saving data for xmgrace ===\n');
baseName = 'ZrNiSn_IIS_900K_scattering';

% Create output directory
outputDir = 'IIS_xmgrace_900K_scattering';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% ===== 1. SAVE PER-BAND INDIVIDUAL K-POINTS =====
fprintf('\n--- Individual k-point files per band ---\n');
for b = 1:numBands
    if h_created(b) && ~isempty(all_E_by_band{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_individual.dat', baseName, b));
        band_data = [all_E_by_band{b}, all_sca_by_band{b}];
        dlmwrite(filename, band_data, 'delimiter', '\t', 'precision', '%.6e');
        fprintf('✓ Band %d: %d individual points → %s\n', b, length(all_E_by_band{b}), filename);
    end
end

% ===== 2. SAVE PER-BAND AVERAGES =====
fprintf('\n--- Per-band average files ---\n');
for b = 1:numBands
    if h_created(b) && ~isempty(E_avg_by_band{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_average.dat', baseName, b));
        
        % Create formatted file with header
        fid = fopen(filename, 'w');
        fprintf(fid, '# ZrNiSn IIS - Band %d\n', b);
        fprintf(fid, '# EF=0, T=900K\n');
        fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
        for j = 1:length(E_avg_by_band{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band{b}(j), sca_avg_by_band{b}(j));
        end
        fclose(fid);
        
        fprintf('✓ Band %d: %d averaged points → %s\n', b, length(E_avg_by_band{b}), filename);
    end
end

% ===== 3. SAVE OVERALL AVERAGE (ACROSS ALL BANDS) =====
fprintf('\n--- Overall average file ---\n');
if ~isempty(E_avg_all)
    filename = fullfile(outputDir, sprintf('%s_overall_average.dat', baseName));
    
    fid = fopen(filename, 'w');
    fprintf(fid, '# ZrNiSn IIS - Overall Average (All Bands)\n');
    fprintf(fid, '# EF=0, T=900K\n');
    fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
    for j = 1:length(E_avg_all)
        fprintf(fid, '%.6f  %.6e\n', E_avg_all(j), sca_avg_all(j));
    end
    fclose(fid);
    
    fprintf('✓ Overall average: %d points → %s\n', length(E_avg_all), filename);
end

% ===== 4. SAVE COMBINED FILE (OPTIONAL) =====
fprintf('\n--- Combined files (optional) ---\n');
if ~isempty(all_E)
    % Combined individual points with band labels
    filename = fullfile(outputDir, sprintf('%s_all_bands_combined.dat', baseName));
    combined_data = [all_E, all_sca, all_band];
    dlmwrite(filename, combined_data, 'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Combined file with band labels: %s\n', filename);
    
    % Combined individual points without band labels
    filename = fullfile(outputDir, sprintf('%s_individual_all.dat', baseName));
    dlmwrite(filename, [all_E, all_sca], 'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Combined file without labels: %s\n', filename);
end

% ===== 5. CREATE XMGRACE BATCH FILE FOR ALL BANDS =====
fprintf('\n--- xmgrace batch file ---\n');
filename = fullfile(outputDir, sprintf('%s_plot.bfile', baseName));
fid = fopen(filename, 'w');

fprintf(fid, '# xmgrace batch file for IIS scattering rates\n');
fprintf(fid, '# Generated from MATLAB\n\n');
fprintf(fid, 'title "ZrNiSn IIS - EF=0, T=900K"\n');
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [fs\\S-1\\N]"\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'legend on\n\n');

% Add per-band averages as separate datasets
set_count = 0;
for b = 1:numBands
    if h_created(b)
        fprintf(fid, 'READ NXY "%s_band%02d_average.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d legend "Band %d"\n', set_count, b);
        set_count = set_count + 1;
    end
end

% Add overall average
fprintf(fid, 'READ NXY "%s_overall_average.dat"\n', baseName);
fprintf(fid, 's%d line color 1\n', set_count);
fprintf(fid, 's%d line linewidth 3\n', set_count);
fprintf(fid, 's%d legend "Overall Average"\n', set_count);

fprintf(fid, '\nAUTOSCALE\n');
fclose(fid);
fprintf('✓ xmgrace batch file → %s\n', filename);

% ===== 6. CREATE SEPARATE BATCH FILES FOR EACH BAND =====
fprintf('\n--- Individual band batch files ---\n');
for b = 1:numBands
    if h_created(b)
        filename = fullfile(outputDir, sprintf('%s_band%02d_plot.bfile', baseName, b));
        fid = fopen(filename, 'w');
        
        fprintf(fid, '# xmgrace batch file for IIS Band %d only\n', b);
        fprintf(fid, 'title "ZrNiSn IIS - Band %d, EF=0, T=900K"\n', b);
        fprintf(fid, 'xaxis label "Energy [eV]"\n');
        fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [fs\\S-1\\N]"\n');
        fprintf(fid, 'yaxis scale logarithmic\n');
        fprintf(fid, 'legend on\n\n');
        
        % Individual k-points (as scatter)
        fprintf(fid, 'READ NXY "%s_band%02d_individual.dat"\n', baseName, b);
        fprintf(fid, 's0 symbol 1\n');  % Circle
        fprintf(fid, 's0 symbol size 0.3\n');
        fprintf(fid, 's0 symbol fill pattern 1\n');
        fprintf(fid, 's0 line type 0\n');  % No line
        fprintf(fid, 's0 legend "Band %d - Individual"\n', b);
        
        % Band average
        fprintf(fid, 'READ NXY "%s_band%02d_average.dat"\n', baseName, b);
        fprintf(fid, 's1 line color %d\n', b);
        fprintf(fid, 's1 line linewidth 2.5\n');
        fprintf(fid, 's1 legend "Band %d - Average"\n', b);
        
        % Overall average (for comparison)
        fprintf(fid, 'READ NXY "%s_overall_average.dat"\n', baseName);
        fprintf(fid, 's2 line color 1\n');
        fprintf(fid, 's2 line linewidth 2\n');
        fprintf(fid, 's2 line style 2\n');  % Dashed
        fprintf(fid, 's2 legend "Overall Average"\n');
        
        fprintf(fid, '\nAUTOSCALE\n');
        fclose(fid);
        fprintf('✓ Band %d batch file → %s\n', b, filename);
    end
end

% ===== 7. CREATE README FILE =====
filename = fullfile(outputDir, 'README.txt');
fid = fopen(filename, 'w');
fprintf(fid, '========================================\n');
fprintf(fid, 'ZrNiSn IIS Scattering Data\n');
fprintf(fid, '========================================\n');
fprintf(fid, 'Parameters: EF=0, T=900K\n');
fprintf(fid, 'Number of bands with data: %d\n', sum(h_created));
fprintf(fid, '\n--- FILE NAMING CONVENTION ---\n');
fprintf(fid, '*_bandXX_individual.dat  : All k-point scattering rates for band XX\n');
fprintf(fid, '*_bandXX_average.dat     : Energy-averaged rates for band XX\n');
fprintf(fid, '*_overall_average.dat    : Average across all bands\n');
fprintf(fid, '*_all_bands_combined.dat : All data with band index column\n');
fprintf(fid, '*_plot.bfile             : xmgrace batch file for all bands\n');
fprintf(fid, '*_bandXX_plot.bfile      : xmgrace batch file for band XX only\n');
fprintf(fid, '\n--- BAND STATISTICS ---\n');
for b = 1:numBands
    if h_created(b)
        fprintf(fid, 'Band %d: %d individual points, %d averages\n', ...
                b, length(all_sca_by_band{b}), length(sca_avg_by_band{b}));
    else
        fprintf(fid, 'Band %d: No data\n', b);
    end
end
fprintf(fid, '\n--- USAGE ---\n');
fprintf(fid, 'To plot all bands:\n');
fprintf(fid, '  xmgrace -batch %s_plot.bfile\n', baseName);
fprintf(fid, '\nTo plot single band (e.g., band 1):\n');
fprintf(fid, '  xmgrace -batch %s_band01_plot.bfile\n', baseName);
fclose(fid);
fprintf('\n✓ README file → %s\n', filename);

%% Summary
fprintf('\n========================================\n');
fprintf('PROCESSING COMPLETE\n');
fprintf('========================================\n');
fprintf('Scattering type: IIS\n');
fprintf('Parameters: EF=0, T=900K\n');
fprintf('Energy points processed: %d\n', max_idx);
fprintf('Bands with data: %d/%d\n', sum(h_created), numBands);
fprintf('\nFiles created in: %s\n', outputDir);
fprintf('  - %d band individual files\n', sum(h_created));
fprintf('  - %d band average files\n', sum(h_created));
fprintf('  - 1 overall average file\n');
fprintf('  - 2 combined files\n');
fprintf('  - 1 all-bands batch file\n');
fprintf('  - %d single-band batch files\n', sum(h_created));
fprintf('  - 1 README file\n');
fprintf('\nTotal files: %d\n', 3*sum(h_created) + 5);
fprintf('========================================\n');
