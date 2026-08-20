%% Load all data files
load 'TE_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_matthiessen_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_screenedPOP_ZrNiSn_kScan_holes.mat';

%% Initialize for plotting
figure;
fig = gcf;
fig.Position(3:4) = [800, 600];
hold on;

numBands = 3;  % Using 3 bands for IIS and POP, 8 for Matthiessen
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
h_created = false(numBands, 1);

% Initialize arrays for data collection
E_avg_all = [];
sca_avg_all = [];

% Initialize per-band storage
E_avg_by_band = cell(numBands, 1);
sca_avg_by_band = cell(numBands, 1);
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);

% For combined storage
all_E = [];
all_sca = [];
all_band = [];

% Configuration for each scattering type
EF_index = 8;
T_index_matt = 1;         % Matthiessen: 1=300K
T_index_IIS = 1;          % IIS: 7=900K
T_index_POP = 1;          % POP: 1=300K

fprintf('Processing TOTAL scattering rate (Matthiessen + IIS + Screened POP)\n');
fprintf('EF = 0 for all mechanisms\n');

%% Main processing loop
max_idx = 451;  % Fixed 451 energy points (used by IIS and POP)

for i = 1:max_idx
    E_row = [];
    sca_row = [];
    
    for b = 1:numBands
        % Initialize total scattering rate for this band/energy
        sca_k_total = [];
        E = state_ID(i,b).E;
        
        % ===== MATTHIESSEN CONTRIBUTION =====
        if b <= 8  % Matthiessen has 8 bands
            if i <= size(state_ID, 1) && ~isempty(taus_matth(i).x)
                % Extract tau for Matthiessen
                if ndims(taus_matth(i).x) == 3
                    tau_k_matt = squeeze(taus_matth(i).x(8, T_index_matt, :));
                elseif ndims(taus_matth(i).x) == 2
                    if size(taus_matth(i).x, 2) >= T_index_matt
                        tau_k_matt = taus_matth(i).x(:, T_index_matt);
                    else
                        tau_k_matt = taus_matth(i).x(:, 1);
                    end
                else
                    tau_k_matt = taus_matth(i).x(:);
                end
                
                positive_mask_matt = tau_k_matt > 0;
                if any(positive_mask_matt)
                    tau_k_matt_positive = tau_k_matt(positive_mask_matt);
                    sca_k_matt = 1 ./ tau_k_matt_positive;
                    
                    % Initialize or add to total
                    if isempty(sca_k_total)
                        sca_k_total = sca_k_matt;
                    else
                        % Interpolate or match k-points (use minimum length)
                        min_len = min(length(sca_k_total), length(sca_k_matt));
                        sca_k_total = sca_k_total(1:min_len) + sca_k_matt(1:min_len);
                    end
                end
            end
        end
        
        % ===== IIS CONTRIBUTION =====
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            tau_3D_IIS = taus_IIS(i,b).x;
            
            if size(tau_3D_IIS, 1) >= EF_index && size(tau_3D_IIS, 2) >= T_index_IIS
                tau_k_IIS = squeeze(tau_3D_IIS(EF_index, T_index_IIS, :));
                positive_mask_IIS = tau_k_IIS > 0;
                
                if any(positive_mask_IIS)
                    tau_k_IIS_positive = tau_k_IIS(positive_mask_IIS);
                    sca_k_IIS = 1 ./ tau_k_IIS_positive;
                    
                    if isempty(sca_k_total)
                        sca_k_total = sca_k_IIS;
                    else
                        min_len = min(length(sca_k_total), length(sca_k_IIS));
                        sca_k_total = sca_k_total(1:min_len) + sca_k_IIS(1:min_len);
                    end
                end
            end
        end
        
        % ===== SCREENED POP CONTRIBUTION =====
        if ~isempty(taus_POP(i,b).x) && size(taus_POP(i,b).x, 3) >= 1
            tau_3D_POP = taus_POP(i,b).x;
            
            if size(tau_3D_POP, 1) >= EF_index && size(tau_3D_POP, 2) >= T_index_POP
                tau_k_POP = squeeze(tau_3D_POP(EF_index, T_index_POP, :));
                positive_mask_POP = tau_k_POP > 0;
                
                if any(positive_mask_POP)
                    tau_k_POP_positive = tau_k_POP(positive_mask_POP);
                    sca_k_POP = 1 ./ tau_k_POP_positive;
                    
                    if isempty(sca_k_total)
                        sca_k_total = sca_k_POP;
                    else
                        min_len = min(length(sca_k_total), length(sca_k_POP));
                        sca_k_total = sca_k_total(1:min_len) + sca_k_POP(1:min_len);
                    end
                end
            end
        end
        
        % ===== PROCESS TOTAL SCATTERING RATE =====
        if ~isempty(sca_k_total) && any(sca_k_total > 0)
            positive_mask_total = sca_k_total > 0;
            sca_k = sca_k_total(positive_mask_total);
            
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
title('ZrNiSn Total Scattering Rate (Matt + IIS + POP)', 'FontSize', 14);

%% Create legend
legend_handles = [];
legend_labels = {};

for b = 1:numBands
    if h_created(b)
        legend_handles = [legend_handles; h(b)];
        legend_labels{end+1} = sprintf('Band %d', b);
    end
end

legend_handles = [legend_handles; avg_line];
legend_labels{end+1} = 'Average';

if ~isempty(legend_handles)
    lgd = legend(legend_handles, legend_labels);
    lgd.FontSize = 10;
    lgd.Location = 'best';
else
    warning('No valid data to create legend');
end

%% ===== SAVE DATA FOR XMGRACE =====
fprintf('\n=== Saving total scattering data for xmgrace ===\n');
baseName = 'ZrNiSn_TotalScattering_Combined';

% Create output directory
outputDir = 'TotalScattering_xmgrace_Combined';
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
        
        fid = fopen(filename, 'w');
        fprintf(fid, '# ZrNiSn Total Scattering - Band %d\n', b);
        fprintf(fid, '# Sum of: Matthiessen(300K) + IIS(900K) + Screened POP(300K)\n');
        fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
        for j = 1:length(E_avg_by_band{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band{b}(j), sca_avg_by_band{b}(j));
        end
        fclose(fid);
        
        fprintf('✓ Band %d: %d averaged points → %s\n', b, length(E_avg_by_band{b}), filename);
    end
end

% ===== 3. SAVE OVERALL AVERAGE =====
fprintf('\n--- Overall average file ---\n');
if ~isempty(E_avg_all)
    filename = fullfile(outputDir, sprintf('%s_overall_average.dat', baseName));
    
    fid = fopen(filename, 'w');
    fprintf(fid, '# ZrNiSn Total Scattering - Overall Average (All Bands)\n');
    fprintf(fid, '# Sum of: Matthiessen(300K) + IIS(900K) + Screened POP(300K)\n');
    fprintf(fid, '# Energy (eV)  Scattering Rate (fs^-1)\n');
    for j = 1:length(E_avg_all)
        fprintf(fid, '%.6f  %.6e\n', E_avg_all(j), sca_avg_all(j));
    end
    fclose(fid);
    
    fprintf('✓ Overall average: %d points → %s\n', length(E_avg_all), filename);
end

% ===== 4. SAVE COMBINED FILE =====
fprintf('\n--- Combined files ---\n');
if ~isempty(all_E)
    % Combined with band labels
    filename = fullfile(outputDir, sprintf('%s_all_bands_combined.dat', baseName));
    combined_data = [all_E, all_sca, all_band];
    dlmwrite(filename, combined_data, 'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Combined file with band labels: %s\n', filename);
    
    % Combined without band labels
    filename = fullfile(outputDir, sprintf('%s_individual_all.dat', baseName));
    dlmwrite(filename, [all_E, all_sca], 'delimiter', '\t', 'precision', '%.6e');
    fprintf('✓ Combined file without labels: %s\n', filename);
end

% ===== 5. CREATE XMGRACE BATCH FILE =====
fprintf('\n--- xmgrace batch file ---\n');
filename = fullfile(outputDir, sprintf('%s_plot.bfile', baseName));
fid = fopen(filename, 'w');

fprintf(fid, '# xmgrace batch file for Total scattering rates\n');
fprintf(fid, '# Generated from MATLAB\n\n');
fprintf(fid, 'title "ZrNiSn Total Scattering (Matt + IIS + POP)"\n');
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [fs\\S-1\\N]"\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'legend on\n\n');

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
        
        fprintf(fid, '# xmgrace batch file for Total Scattering Band %d only\n', b);
        fprintf(fid, 'title "ZrNiSn Total Scattering - Band %d"\n', b);
        fprintf(fid, 'xaxis label "Energy [eV]"\n');
        fprintf(fid, 'yaxis label "Scattering Rate \\xt\\f{Symbol}\\f{Times-Italic}\\f{} \\S-1\\N [fs\\S-1\\N]"\n');
        fprintf(fid, 'yaxis scale logarithmic\n');
        fprintf(fid, 'legend on\n\n');
        
        % Individual k-points
        fprintf(fid, 'READ NXY "%s_band%02d_individual.dat"\n', baseName, b);
        fprintf(fid, 's0 symbol 1\n');
        fprintf(fid, 's0 symbol size 0.3\n');
        fprintf(fid, 's0 symbol fill pattern 1\n');
        fprintf(fid, 's0 line type 0\n');
        fprintf(fid, 's0 legend "Band %d - Individual"\n', b);
        
        % Band average
        fprintf(fid, 'READ NXY "%s_band%02d_average.dat"\n', baseName, b);
        fprintf(fid, 's1 line color %d\n', b);
        fprintf(fid, 's1 line linewidth 2.5\n');
        fprintf(fid, 's1 legend "Band %d - Average"\n', b);
        
        % Overall average
        fprintf(fid, 'READ NXY "%s_overall_average.dat"\n', baseName);
        fprintf(fid, 's2 line color 1\n');
        fprintf(fid, 's2 line linewidth 2\n');
        fprintf(fid, 's2 line style 2\n');
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
fprintf(fid, 'ZrNiSn Total Scattering Rate Data\n');
fprintf(fid, '========================================\n');
fprintf(fid, 'Scattering mechanisms included:\n');
fprintf(fid, '  1. Matthiessen (T = 300K)\n');
fprintf(fid, '  2. Ionized Impurity Scattering - IIS (T = 900K)\n');
fprintf(fid, '  3. Screened Polar Optical Phonon - POP (T = 300K)\n');
fprintf(fid, '\nFormula: tau_total^-1 = tau_matt^-1 + tau_IIS^-1 + tau_POP^-1\n');
fprintf(fid, 'EF = 0 eV for all mechanisms\n');
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
fprintf('Total Scattering Rate = Matt + IIS + POP\n');
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
