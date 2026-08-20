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
% Simple (unweighted) averages
E_avg_all_simple = [];
sca_avg_all_simple = [];

% DOS-weighted averages
E_avg_all_DOS = [];
sca_avg_all_DOS = [];

% Initialize per-band storage
% Simple averages
E_avg_by_band_simple = cell(numBands, 1);
sca_avg_by_band_simple = cell(numBands, 1);

% DOS-weighted averages
E_avg_by_band_DOS = cell(numBands, 1);
sca_avg_by_band_DOS = cell(numBands, 1);

% Individual k-points
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);
all_DOS_by_band = cell(numBands, 1);  % Store DOS values

% For combined storage
all_E = [];
all_sca = [];
all_band = [];
all_DOS = [];

% Configuration
EF_index = 8;            % Index for EF=0
T_index = 1;             % Temperature index (1=300K, adjust if needed)

fprintf('Processing IIS data for EF=0, T=300K (BOTH simple and DOS-weighted)\n');
fprintf('Simple mean: S_simple = mean(S_k)\n');
fprintf('DOS-weighted: S_DOS = sum(S_k * DOS_k) / sum(DOS_k)\n');

%% Main processing loop
max_idx = 451;  % Fixed 451 energy points

for i = 1:max_idx
    % Arrays for this energy point (across bands)
    E_row = [];
    sca_simple_row = [];  % Simple mean per band
    sca_DOS_row = [];     % DOS-weighted mean per band
    DOS_row = [];         % Total DOS per band
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            % Get energy for this point
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= EF_index && size(tau_3D, 2) >= T_index
                % Extract tau(EF=0, T=300K, all k-points)
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                
                % Extract DOS for all k-points
                if isfield(state_ID(i,b), 'DOS') && ~isempty(state_ID(i,b).DOS)
                    DOS_k = state_ID(i,b).DOS(:);  % Ensure column vector
                    
                    % Make sure DOS_k has same length as tau_k
                    if length(DOS_k) ~= length(tau_k)
                        warning('DOS and tau size mismatch at E=%.4f, band %d. Using uniform weights.', E, b);
                        DOS_k = ones(size(tau_k));  % Fall back to uniform weights
                    end
                else
                    warning('No DOS data at E=%.4f, band %d. Using uniform weights.', E, b);
                    DOS_k = ones(size(tau_k));  % Uniform weights if no DOS
                end
                
                % Filter: Keep only positive tau values
                positive_mask = tau_k > 0;
                
                if any(positive_mask)
                    tau_k_positive = tau_k(positive_mask);
                    DOS_k_positive = DOS_k(positive_mask);
                    sca_k = 1 ./ tau_k_positive;  % Scattering rates
                    
                    % Store individual k-points for this band (unweighted)
                    n_points = length(sca_k);
                    all_E_by_band{b} = [all_E_by_band{b}; repmat(E, n_points, 1)];
                    all_sca_by_band{b} = [all_sca_by_band{b}; sca_k];
                    all_DOS_by_band{b} = [all_DOS_by_band{b}; DOS_k_positive];
                    
                    % Also store in combined arrays
                    all_E = [all_E; repmat(E, n_points, 1)];
                    all_sca = [all_sca; sca_k];
                    all_band = [all_band; repmat(b, n_points, 1)];
                    all_DOS = [all_DOS; DOS_k_positive];
                    
                    % Scatter plot for all k-points at this energy
                    htmp = scatter(repmat(E, size(sca_k)), sca_k, 30, palette(b,:), 'filled');
                    
                    % SIMPLE mean over k-points for this band
                    sca_band_avg_simple = mean(sca_k);
                    
                    % DOS-WEIGHTED mean over k-points for this band
                    sca_band_avg_DOS = sum(sca_k .* DOS_k_positive) / sum(DOS_k_positive);
                    
                    % Store handle only if this is the first valid data for this band
                    if ~h_created(b)
                        h(b) = htmp;
                        h_created(b) = true;
                    end
                    
                    % Plot both averages for this energy
                    % Simple mean: filled circle
                    scatter(E, sca_band_avg_simple, 60, palette(b,:), 'o', 'LineWidth', 1.5);
                    
                    % DOS-weighted mean: filled diamond
                    scatter(E, sca_band_avg_DOS, 60, palette(b,:), 'd', 'LineWidth', 1.5);
                    
                    % Store both types of averages for this band
                    E_avg_by_band_simple{b}(end+1) = E;
                    sca_avg_by_band_simple{b}(end+1) = sca_band_avg_simple;
                    
                    E_avg_by_band_DOS{b}(end+1) = E;
                    sca_avg_by_band_DOS{b}(end+1) = sca_band_avg_DOS;
                    
                    % Collect for total averages over bands
                    E_row(end+1) = E;
                    sca_simple_row(end+1) = sca_band_avg_simple;
                    sca_DOS_row(end+1) = sca_band_avg_DOS;
                    DOS_row(end+1) = sum(DOS_k_positive);  % Total DOS for this band
                end
            end
        end
    end
    
    % Calculate overall averages across bands
    if ~isempty(E_row)
        % Simple average across bands
        E_avg_all_simple(end+1) = mean(E_row);
        sca_avg_all_simple(end+1) = mean(sca_simple_row);
        
        % DOS-weighted average across bands
        if ~isempty(DOS_row) && sum(DOS_row) > 0
            E_avg_all_DOS(end+1) = sum(E_row .* DOS_row) / sum(DOS_row);
            sca_avg_all_DOS(end+1) = sum(sca_DOS_row .* DOS_row) / sum(DOS_row);
        else
            E_avg_all_DOS(end+1) = mean(E_row);
            sca_avg_all_DOS(end+1) = mean(sca_DOS_row);
        end
    end
    
    % Progress indicator
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', i, max_idx);
    end
end

%% Overlay average lines
% Simple average line (black solid)
avg_line_simple = plot(E_avg_all_simple, sca_avg_all_simple, 'k-', 'LineWidth', 2.5);

% DOS-weighted average line (red dashed)
avg_line_DOS = plot(E_avg_all_DOS, sca_avg_all_DOS, 'r--', 'LineWidth', 2.5);

%% Format the plot
set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([0, 0.45]);
grid on;

xlabel('Energy [eV]', 'FontSize', 12);
ylabel('\tau^{-1} [fs^{-1}]', 'FontSize', 12);
title('ZrNiSn IIS - EF=0, T=300K (Simple vs DOS-weighted)', 'FontSize', 14);

%% Add annotation for symbol meaning
text(0.02, 0.98, '○ = Simple mean per band', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
text(0.02, 0.94, '◇ = DOS-weighted mean per band', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
text(0.02, 0.90, '— = Overall simple average', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
text(0.02, 0.86, '- - = Overall DOS-weighted average', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);

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

% Add both average lines
legend_handles = [legend_handles; avg_line_simple; avg_line_DOS];
legend_labels{end+1} = 'Overall Simple Avg';
legend_labels{end+1} = 'Overall DOS-Weighted Avg';

% Only create legend if we have valid handles
if ~isempty(legend_handles)
    lgd = legend(legend_handles, legend_labels);
    lgd.FontSize = 9;
    lgd.Location = 'best';
else
    warning('No valid data to create legend');
end

%% ===== SAVE DATA FOR XMGRACE - BOTH SIMPLE AND DOS-WEIGHTED =====
fprintf('\n=== Saving BOTH simple and DOS-weighted data ===\n');
baseName = 'ZrNiSn_IIS_300K';

% Create output directory
outputDir = 'IIS_xmgrace_300K_both';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% ===== 1. SAVE PER-BAND INDIVIDUAL K-POINTS (with DOS) =====
fprintf('\n--- Individual k-point files per band (with DOS) ---\n');
for b = 1:numBands
    if h_created(b) && ~isempty(all_E_by_band{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_individual.dat', baseName, b));
        band_data = [all_E_by_band{b}, all_sca_by_band{b}, all_DOS_by_band{b}];
        dlmwrite(filename, band_data, 'delimiter', '\t', 'precision', '%.6e');
        fprintf('✓ Band %d: %d individual points (with DOS) → %s\n', b, length(all_E_by_band{b}), filename);
    end
end

% ===== 2. SAVE PER-BAND SIMPLE AVERAGES =====
fprintf('\n--- Per-band SIMPLE average files ---\n');
for b = 1:numBands
    if h_created(b) && ~isempty(E_avg_by_band_simple{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_simple_average.dat', baseName, b));
        
        fid = fopen(filename, 'w');
        fprintf(fid, '# ZrNiSn IIS - Band %d (SIMPLE mean)\n', b);
        fprintf(fid, '# EF=0, T=300K\n');
        fprintf(fid, '# Weighting: S_simple = mean(S_k)\n');
        fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
        for j = 1:length(E_avg_by_band_simple{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_simple{b}(j), sca_avg_by_band_simple{b}(j));
        end
        fclose(fid);
        
        fprintf('✓ Band %d: %d simple averaged points → %s\n', b, length(E_avg_by_band_simple{b}), filename);
    end
end

% ===== 3. SAVE PER-BAND DOS-WEIGHTED AVERAGES =====
fprintf('\n--- Per-band DOS-WEIGHTED average files ---\n');
for b = 1:numBands
    if h_created(b) && ~isempty(E_avg_by_band_DOS{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_DOSweighted_average.dat', baseName, b));
        
        fid = fopen(filename, 'w');
        fprintf(fid, '# ZrNiSn IIS - Band %d (DOS-WEIGHTED mean)\n', b);
        fprintf(fid, '# EF=0, T=300K\n');
        fprintf(fid, '# Weighting: S_DOS = sum(S_k * DOS_k) / sum(DOS_k)\n');
        fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
        for j = 1:length(E_avg_by_band_DOS{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_DOS{b}(j), sca_avg_by_band_DOS{b}(j));
        end
        fclose(fid);
        
        fprintf('✓ Band %d: %d DOS-weighted averaged points → %s\n', b, length(E_avg_by_band_DOS{b}), filename);
    end
end

% ===== 4. SAVE OVERALL SIMPLE AVERAGE =====
fprintf('\n--- Overall SIMPLE average file ---\n');
if ~isempty(E_avg_all_simple)
    filename = fullfile(outputDir, sprintf('%s_overall_simple_average.dat', baseName));
    
    fid = fopen(filename, 'w');
    fprintf(fid, '# ZrNiSn IIS - Overall SIMPLE Average (All Bands)\n');
    fprintf(fid, '# EF=0, T=300K\n');
    fprintf(fid, '# Weighting: Simple mean across bands\n');
    fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
    for j = 1:length(E_avg_all_simple)
        fprintf(fid, '%.6f  %.6e\n', E_avg_all_simple(j), sca_avg_all_simple(j));
    end
    fclose(fid);
    
    fprintf('✓ Overall simple average: %d points → %s\n', length(E_avg_all_simple), filename);
end

% ===== 5. SAVE OVERALL DOS-WEIGHTED AVERAGE =====
fprintf('\n--- Overall DOS-WEIGHTED average file ---\n');
if ~isempty(E_avg_all_DOS)
    filename = fullfile(outputDir, sprintf('%s_overall_DOSweighted_average.dat', baseName));
    
    fid = fopen(filename, 'w');
    fprintf(fid, '# ZrNiSn IIS - Overall DOS-WEIGHTED Average (All Bands)\n');
    fprintf(fid, '# EF=0, T=300K\n');
    fprintf(fid, '# Weighting: S_total = sum(S_band * DOS_band) / sum(DOS_band)\n');
    fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
    for j = 1:length(E_avg_all_DOS)
        fprintf(fid, '%.6f  %.6e\n', E_avg_all_DOS(j), sca_avg_all_DOS(j));
    end
    fclose(fid);
    
    fprintf('✓ Overall DOS-weighted average: %d points → %s\n', length(E_avg_all_DOS), filename);
end

% ===== 6. SAVE COMPARISON FILE (BOTH AVERAGES) =====
fprintf('\n--- Comparison file (simple vs DOS-weighted) ---\n');
if ~isempty(E_avg_all_simple) && ~isempty(E_avg_all_DOS)
    filename = fullfile(outputDir, sprintf('%s_comparison.dat', baseName));
    
    fid = fopen(filename, 'w');
    fprintf(fid, '# ZrNiSn IIS - Comparison of Simple vs DOS-Weighted\n');
    fprintf(fid, '# EF=0, T=300K\n');
    fprintf(fid, '# Energy  Simple_Avg  DOS-Weighted_Avg  Difference  Ratio(DOS/Simple)\n');
    for j = 1:min(length(E_avg_all_simple), length(E_avg_all_DOS))
        diff = sca_avg_all_DOS(j) - sca_avg_all_simple(j);
        ratio = sca_avg_all_DOS(j) / sca_avg_all_simple(j);
        fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.4f\n', ...
                E_avg_all_simple(j), sca_avg_all_simple(j), sca_avg_all_DOS(j), diff, ratio);
    end
    fclose(fid);
    
    fprintf('✓ Comparison file: %d points → %s\n', j, filename);
end

% ===== 7. SAVE COMBINED FILES WITH DOS INFORMATION =====
fprintf('\n--- Combined files with DOS ---\n');
if ~isempty(all_E)
    % Combined individual points with band labels and DOS
    filename = fullfile(outputDir, sprintf('%s_all_bands_combined.dat', baseName));
    combined_data = [all_E, all_sca, all_DOS, all_band];
    dlmwrite(filename, combined_data, 'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Combined file with band labels and DOS: %s\n', filename);
    
    % Combined individual points without band labels (but with DOS)
    filename = fullfile(outputDir, sprintf('%s_individual_all.dat', baseName));
    dlmwrite(filename, [all_E, all_sca, all_DOS], 'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Combined file with DOS (no labels): %s\n', filename);
end

% ===== 8. CREATE XMGRACE BATCH FILE FOR ALL BANDS =====
fprintf('\n--- xmgrace batch file ---\n');
filename = fullfile(outputDir, sprintf('%s_plot.bfile', baseName));
fid = fopen(filename, 'w');

fprintf(fid, '# xmgrace batch file for IIS scattering rates (Simple vs DOS-weighted)\n');
fprintf(fid, '# Generated from MATLAB\n\n');
fprintf(fid, 'title "ZrNiSn IIS - EF=0, T=300K (Simple vs DOS-Weighted)"\n');
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [fs\\S-1\\N]"\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'legend on\n\n');

% Add per-band simple averages
set_count = 0;
for b = 1:numBands
    if h_created(b)
        fprintf(fid, 'READ NXY "%s_band%02d_simple_average.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d line style 1\n', set_count);  % Solid
        fprintf(fid, 's%d legend "Band %d (Simple)"\n', set_count, b);
        set_count = set_count + 1;
    end
end

% Add per-band DOS-weighted averages
for b = 1:numBands
    if h_created(b)
        fprintf(fid, 'READ NXY "%s_band%02d_DOSweighted_average.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d line style 2\n', set_count);  % Dashed
        fprintf(fid, 's%d legend "Band %d (DOS-wt)"\n', set_count, b);
        set_count = set_count + 1;
    end
end

% Add overall simple average
fprintf(fid, 'READ NXY "%s_overall_simple_average.dat"\n', baseName);
fprintf(fid, 's%d line color 1\n', set_count);
fprintf(fid, 's%d line linewidth 3\n', set_count);
fprintf(fid, 's%d line style 1\n', set_count);  % Solid
fprintf(fid, 's%d legend "Overall Simple Avg"\n', set_count);
set_count = set_count + 1;

% Add overall DOS-weighted average
fprintf(fid, 'READ NXY "%s_overall_DOSweighted_average.dat"\n', baseName);
fprintf(fid, 's%d line color 2\n', set_count);
fprintf(fid, 's%d line linewidth 3\n', set_count);
fprintf(fid, 's%d line style 2\n', set_count);  % Dashed
fprintf(fid, 's%d legend "Overall DOS-Weighted Avg"\n', set_count);

fprintf(fid, '\nAUTOSCALE\n');
fclose(fid);
fprintf('✓ xmgrace batch file → %s\n', filename);

% ===== 9. CREATE SEPARATE BATCH FILES FOR EACH BAND =====
fprintf('\n--- Individual band batch files ---\n');
for b = 1:numBands
    if h_created(b)
        filename = fullfile(outputDir, sprintf('%s_band%02d_plot.bfile', baseName, b));
        fid = fopen(filename, 'w');
        
        fprintf(fid, '# xmgrace batch file for IIS Band %d (Simple vs DOS-weighted)\n', b);
        fprintf(fid, 'title "ZrNiSn IIS - Band %d, EF=0, T=300K"\n', b);
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
        
        % Simple average
        fprintf(fid, 'READ NXY "%s_band%02d_simple_average.dat"\n', baseName, b);
        fprintf(fid, 's1 line color %d\n', b);
        fprintf(fid, 's1 line linewidth 2.5\n');
        fprintf(fid, 's1 line style 1\n');  % Solid
        fprintf(fid, 's1 legend "Band %d - Simple Avg"\n', b);
        
        % DOS-weighted average
        fprintf(fid, 'READ NXY "%s_band%02d_DOSweighted_average.dat"\n', baseName, b);
        fprintf(fid, 's2 line color %d\n', b);
        fprintf(fid, 's2 line linewidth 2.5\n');
        fprintf(fid, 's2 line style 2\n');  % Dashed
        fprintf(fid, 's2 legend "Band %d - DOS-Weighted Avg"\n', b);
        
        % Overall averages for comparison
        fprintf(fid, 'READ NXY "%s_overall_simple_average.dat"\n', baseName);
        fprintf(fid, 's3 line color 1\n');
        fprintf(fid, 's3 line linewidth 2\n');
        fprintf(fid, 's3 line style 1\n');
        fprintf(fid, 's3 legend "Overall Simple Avg"\n');
        
        fprintf(fid, 'READ NXY "%s_overall_DOSweighted_average.dat"\n', baseName);
        fprintf(fid, 's4 line color 2\n');
        fprintf(fid, 's4 line linewidth 2\n');
        fprintf(fid, 's4 line style 2\n');
        fprintf(fid, 's4 legend "Overall DOS-Weighted Avg"\n');
        
        fprintf(fid, '\nAUTOSCALE\n');
        fclose(fid);
        fprintf('✓ Band %d batch file → %s\n', b, filename);
    end
end

% ===== 10. CREATE README FILE =====
filename = fullfile(outputDir, 'README.txt');
fid = fopen(filename, 'w');
fprintf(fid, '========================================\n');
fprintf(fid, 'ZrNiSn IIS Scattering Data\n');
fprintf(fid, 'BOTH Simple and DOS-Weighted Averages\n');
fprintf(fid, '========================================\n');
fprintf(fid, 'Parameters: EF=0, T=300K\n');
fprintf(fid, 'Number of bands with data: %d\n', sum(h_created));
fprintf(fid, '\n--- AVERAGING METHODS ---\n');
fprintf(fid, '1. SIMPLE mean: S_simple = mean(S_k)\n');
fprintf(fid, '2. DOS-WEIGHTED mean: S_DOS = sum(S_k * DOS_k) / sum(DOS_k)\n');
fprintf(fid, '\n--- FILE NAMING CONVENTION ---\n');
fprintf(fid, '*_bandXX_individual.dat          : All k-point scattering rates + DOS\n');
fprintf(fid, '  Columns: Energy, Scattering Rate, DOS\n');
fprintf(fid, '*_bandXX_simple_average.dat      : SIMPLE averages for band XX\n');
fprintf(fid, '*_bandXX_DOSweighted_average.dat : DOS-WEIGHTED averages for band XX\n');
fprintf(fid, '*_overall_simple_average.dat     : Overall SIMPLE average\n');
fprintf(fid, '*_overall_DOSweighted_average.dat: Overall DOS-WEIGHTED average\n');
fprintf(fid, '*_comparison.dat                 : Comparison of both methods\n');
fprintf(fid, '  Columns: Energy, Simple Avg, DOS-wt Avg, Difference, Ratio\n');
fprintf(fid, '*_all_bands_combined.dat         : All data with DOS and band index\n');
fprintf(fid, '*_plot.bfile                     : xmgrace batch file for all\n');
fprintf(fid, '*_bandXX_plot.bfile              : xmgrace batch file for band XX\n');
fprintf(fid, '\n--- BAND STATISTICS ---\n');
for b = 1:numBands
    if h_created(b)
        fprintf(fid, 'Band %d:\n', b);
        fprintf(fid, '  Individual points: %d\n', length(all_sca_by_band{b}));
        fprintf(fid, '  Simple averages: %d\n', length(sca_avg_by_band_simple{b}));
        fprintf(fid, '  DOS-weighted averages: %d\n', length(sca_avg_by_band_DOS{b}));
        fprintf(fid, '  DOS range: %.2e to %.2e\n', min(all_DOS_by_band{b}), max(all_DOS_by_band{b}));
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
fprintf('Parameters: EF=0, T=300K\n');
fprintf('BOTH simple and DOS-weighted averages calculated\n');
fprintf('Energy points processed: %d\n', max_idx);
fprintf('Bands with data: %d/%d\n', sum(h_created), numBands);
fprintf('\nFiles created in: %s\n', outputDir);
fprintf('  - %d band individual files (with DOS)\n', sum(h_created));
fprintf('  - %d band simple average files\n', sum(h_created));
fprintf('  - %d band DOS-weighted average files\n', sum(h_created));
fprintf('  - 1 overall simple average file\n');
fprintf('  - 1 overall DOS-weighted average file\n');
fprintf('  - 1 comparison file\n');
fprintf('  - 2 combined files\n');
fprintf('  - 1 all-bands batch file\n');
fprintf('  - %d single-band batch files\n', sum(h_created));
fprintf('  - 1 README file\n');
fprintf('\nTotal files: %d\n', 5*sum(h_created) + 7);
fprintf('========================================\n');