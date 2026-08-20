%% Load the data files
%load 'TE_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';

%% Initialize for plotting - Create two figures
% Figure 1: Scattering Rates
figure1 = figure('Name', 'Scattering Rates');
fig1 = gcf;
fig1.Position(3:4) = [1000, 800];
hold on;

% Figure 2: Relaxation Times
figure2 = figure('Name', 'Relaxation Times');
fig2 = gcf;
fig2.Position(3:4) = [1000, 800];
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

% Initialize arrays for data collection
% ===== SCATTERING RATES (τ⁻¹) =====
% Simple averages
E_avg_all_simple_sca = [];
sca_avg_all_simple = [];

% DOS-weighted averages (Method 1: weight S_k directly)
E_avg_all_DOS_sca = [];
sca_avg_all_DOS = [];

% DOS-weighted averages (Method 2: weight τ, then invert)
E_avg_all_DOS_tau_inv = [];
sca_avg_all_DOS_tau_inv = [];

% ===== RELAXATION TIMES (τ) =====
% Simple averages
E_avg_all_simple_tau = [];
tau_avg_all_simple = [];

% DOS-weighted averages (Method 1: weight τ directly)
E_avg_all_DOS_tau = [];
tau_avg_all_DOS = [];

% DOS-weighted averages (Method 2: weight S, then invert)
E_avg_all_DOS_sca_inv = [];
tau_avg_all_DOS_sca_inv = [];

% Initialize per-band storage
% ===== SCATTERING RATES =====
% Simple averages
E_avg_by_band_simple = cell(numBands, 1);
sca_avg_by_band_simple = cell(numBands, 1);

% DOS-weighted (Method 1: weight S)
E_avg_by_band_DOS_sca = cell(numBands, 1);
sca_avg_by_band_DOS = cell(numBands, 1);

% DOS-weighted (Method 2: weight τ, then invert)
E_avg_by_band_DOS_tau_inv = cell(numBands, 1);
sca_avg_by_band_DOS_tau_inv = cell(numBands, 1);

% ===== RELAXATION TIMES =====
% Simple averages
E_avg_by_band_simple_tau = cell(numBands, 1);
tau_avg_by_band_simple = cell(numBands, 1);

% DOS-weighted (Method 1: weight τ)
E_avg_by_band_DOS_tau = cell(numBands, 1);
tau_avg_by_band_DOS = cell(numBands, 1);

% DOS-weighted (Method 2: weight S, then invert)
E_avg_by_band_DOS_sca_inv = cell(numBands, 1);
tau_avg_by_band_DOS_sca_inv = cell(numBands, 1);

% Individual k-points
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);
all_tau_by_band = cell(numBands, 1);
all_DOS_by_band = cell(numBands, 1);

% For combined storage
all_E = [];
all_sca = [];
all_tau = [];
all_band = [];
all_DOS = [];

% Configuration
EF_index = 8;            % Index for EF=0
T_index = 1;             % Temperature index (1=300K, adjust if needed)

fprintf('========================================\n');
fprintf('Processing IIS data for EF=0, T=300K\n');
fprintf('Calculating BOTH scattering rates AND relaxation times\n');
fprintf('========================================\n');
fprintf('\nAveraging Methods:\n');
fprintf('1. Simple mean (unweighted)\n');
fprintf('2. DOS-weighted (Method 1: weight the quantity directly)\n');
fprintf('3. DOS-weighted (Method 2: weight inverse, then invert)\n');

%% Main processing loop
max_idx = 451;  % Fixed 451 energy points

for i = 1:max_idx
    % Arrays for this energy point (across bands)
    E_row = [];
    
    % Scattering rates
    sca_simple_row = [];
    sca_DOS_row = [];
    sca_DOS_tau_inv_row = [];
    
    % Relaxation times
    tau_simple_row = [];
    tau_DOS_row = [];
    tau_DOS_sca_inv_row = [];
    
    % DOS for overall weighting
    DOS_row = [];
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            % Get energy for this point
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= EF_index && size(tau_3D, 2) >= T_index
                % Extract tau(EF=0, T=300K, all k-points)
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                                
                % Extract velocity squared (V²) for all k-points
                if isfield(state_ID(i,b), 'V')
                    % Velocity is stored as [vx, vy, vz] for each k-point
                    V_k = state_ID(i,b).V(:);
		    V2_k = V_k.^2
		else 
                    warning('No velocity data at E=%.4f, band %d. Using uniform V^2.', E, b);
                    V2_k = ones(size(tau_k));
                end

                % Ensure V2_k matches tau_k in size
                if length(V2_k) ~= length(tau_k)
                    warning('Velocity and tau size mismatch at E=%.4f, band %d.', E, b);
                    V2_k = ones(size(tau_k));
                end                

                % Extract DOS for all k-points
                if isfield(state_ID(i,b), 'DOS') && ~isempty(state_ID(i,b).DOS)
                    DOS_k = state_ID(i,b).DOS(:);  % Ensure column vector
                    
                    % Make sure DOS_k has same length as tau_k
                    if length(DOS_k) ~= length(tau_k)
                        warning('DOS and tau size mismatch at E=%.4f, band %d.', E, b);
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
		    V2_k_positive = V2_k(positive_mask);
                    DOS_k_positive = DOS_k(positive_mask);
                    sca_k = 1 ./ tau_k_positive;  % Scattering rates

                    % ==== PHYSICALLY CORRECT: Average scattering rates, then invert ====
                    % Transport weight = V² * g
                    transport_weight = V2_k_positive .* DOS_k_positive;
                    
                    % Store individual k-points for this band
                    n_points = length(sca_k);
                    all_E_by_band{b} = [all_E_by_band{b}; repmat(E, n_points, 1)];
                    all_sca_by_band{b} = [all_sca_by_band{b}; sca_k];
                    all_tau_by_band{b} = [all_tau_by_band{b}; tau_k_positive];
                    all_DOS_by_band{b} = [all_DOS_by_band{b}; DOS_k_positive];
                    
                    % Also store in combined arrays
                    all_E = [all_E; repmat(E, n_points, 1)];
                    all_sca = [all_sca; sca_k];
                    all_tau = [all_tau; tau_k_positive];
                    all_band = [all_band; repmat(b, n_points, 1)];
                    all_DOS = [all_DOS; DOS_k_positive];
                    
                    % ===== SCATTERING RATE AVERAGES =====
                    % 1. Simple mean of scattering rates
                    sca_band_avg_simple = mean(sca_k);
                    
                    % 2. DOS-weighted scattering (Method 1: weight S_k directly)
                    sca_band_avg_DOS = sum(sca_k .* DOS_k_positive) / sum(DOS_k_positive);
                    
                    % 3. DOS-weighted scattering (Method 2: weight τ, then invert)
                    tau_band_avg_DOS = sum(tau_k_positive .* DOS_k_positive) / sum(DOS_k_positive);
                    sca_band_avg_DOS_tau_inv = 1 / tau_band_avg_DOS;
                    
                    % ===== RELAXATION TIME AVERAGES =====
                    % 1. Simple mean of relaxation times
                    tau_band_avg_simple = 1/sca_band_avg_simple;
                    
                    % 2. DOS-weighted relaxation (Method 1: weight τ directly)
                    tau_band_avg_DOS_direct = tau_band_avg_DOS;  % Already calculated above
                    
                    % 3. DOS-weighted relaxation (Method 2: weight S, then invert)
                    sca_band_avg_DOS_direct = sca_band_avg_DOS;  % Already calculated above
                    tau_band_avg_DOS_sca_inv = 1 / sca_band_avg_DOS_direct;
                    
                    % Store scattering rate averages
                    E_avg_by_band_simple{b}(end+1) = E;
                    sca_avg_by_band_simple{b}(end+1) = sca_band_avg_simple;
                    
                    E_avg_by_band_DOS_sca{b}(end+1) = E;
                    sca_avg_by_band_DOS{b}(end+1) = sca_band_avg_DOS;
                    
                    E_avg_by_band_DOS_tau_inv{b}(end+1) = E;
                    sca_avg_by_band_DOS_tau_inv{b}(end+1) = sca_band_avg_DOS_tau_inv;
                    
                    % Store relaxation time averages
                    E_avg_by_band_simple_tau{b}(end+1) = E;
                    tau_avg_by_band_simple{b}(end+1) = tau_band_avg_simple;
                    
                    E_avg_by_band_DOS_tau{b}(end+1) = E;
                    tau_avg_by_band_DOS{b}(end+1) = tau_band_avg_DOS_direct;
                    
                    E_avg_by_band_DOS_sca_inv{b}(end+1) = E;
                    tau_avg_by_band_DOS_sca_inv{b}(end+1) = tau_band_avg_DOS_sca_inv;
                    
                    % Collect for total averages over bands
                    E_row(end+1) = E;
                    sca_simple_row(end+1) = sca_band_avg_simple;
                    sca_DOS_row(end+1) = sca_band_avg_DOS;
                    sca_DOS_tau_inv_row(end+1) = sca_band_avg_DOS_tau_inv;
                    tau_simple_row(end+1) = tau_band_avg_simple;
                    tau_DOS_row(end+1) = tau_band_avg_DOS_direct;
                    tau_DOS_sca_inv_row(end+1) = tau_band_avg_DOS_sca_inv;
                    DOS_row(end+1) = sum(DOS_k_positive);  % Total DOS for this band
                end
            end
        end
    end
    
    % Calculate overall averages across bands
    if ~isempty(E_row)
        % Simple averages across bands (unweighted)
        E_avg_all_simple_sca(end+1) = mean(E_row);
        sca_avg_all_simple(end+1) = mean(sca_simple_row);
        E_avg_all_simple_tau(end+1) = mean(E_row);
        tau_avg_all_simple(end+1) = mean(tau_simple_row);
        
        % DOS-weighted averages across bands
        if ~isempty(DOS_row) && sum(DOS_row) > 0
            E_avg_all_DOS_sca(end+1) = sum(E_row .* DOS_row) / sum(DOS_row);
            sca_avg_all_DOS(end+1) = sum(sca_DOS_row .* DOS_row) / sum(DOS_row);
            sca_avg_all_DOS_tau_inv(end+1) = sum(sca_DOS_tau_inv_row .* DOS_row) / sum(DOS_row);
            
            E_avg_all_DOS_tau(end+1) = sum(E_row .* DOS_row) / sum(DOS_row);
            tau_avg_all_DOS(end+1) = sum(tau_DOS_row .* DOS_row) / sum(DOS_row);
            tau_avg_all_DOS_sca_inv(end+1) = sum(tau_DOS_sca_inv_row .* DOS_row) / sum(DOS_row);
        else
            E_avg_all_DOS_sca(end+1) = mean(E_row);
            sca_avg_all_DOS(end+1) = mean(sca_DOS_row);
            sca_avg_all_DOS_tau_inv(end+1) = mean(sca_DOS_tau_inv_row);
            
            E_avg_all_DOS_tau(end+1) = mean(E_row);
            tau_avg_all_DOS(end+1) = mean(tau_DOS_row);
            tau_avg_all_DOS_sca_inv(end+1) = mean(tau_DOS_sca_inv_row);
        end
    end
    
    % Progress indicator
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', i, max_idx);
    end
end

%% ===== FIGURE 1: SCATTERING RATES =====
figure(figure1);

% Plot individual k-points as small dots for first band only (to avoid clutter)
if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_sca_by_band{1}, 5, [0.7 0.7 0.7], 'filled', ...
            'DisplayName', 'Individual k-points (Band 1)');
end

% Plot per-band averages
h_sca = gobjects(numBands, 3);
legend_entries_sca = {};
legend_count = 0;

for b = 1:numBands
    if ~isempty(E_avg_by_band_simple{b})
        % Simple average (solid line)
        h_sca(b,1) = plot(E_avg_by_band_simple{b}, sca_avg_by_band_simple{b}, '-', ...
                         'Color', palette(b,:), 'LineWidth', 1.5);
        legend_count = legend_count + 1;
        legend_entries_sca{legend_count} = sprintf('Band %d - Simple', b);
        
        % DOS-weighted Method 1 (dashed line)
        h_sca(b,2) = plot(E_avg_by_band_DOS_sca{b}, sca_avg_by_band_DOS{b}, '--', ...
                         'Color', palette(b,:), 'LineWidth', 2);
        legend_count = legend_count + 1;
        legend_entries_sca{legend_count} = sprintf('Band %d - DOS-wt (Method 1)', b);
        
        % DOS-weighted Method 2 (dotted line)
        h_sca(b,3) = plot(E_avg_by_band_DOS_tau_inv{b}, sca_avg_by_band_DOS_tau_inv{b}, ':', ...
                         'Color', palette(b,:), 'LineWidth', 2);
        legend_count = legend_count + 1;
        legend_entries_sca{legend_count} = sprintf('Band %d - DOS-wt (Method 2)', b);
    end
end

% Plot overall averages
if ~isempty(E_avg_all_simple_sca)
    plot(E_avg_all_simple_sca, sca_avg_all_simple, 'k-', 'LineWidth', 3, ...
         'DisplayName', 'Overall Simple');
    plot(E_avg_all_DOS_sca, sca_avg_all_DOS, 'r--', 'LineWidth', 3, ...
         'DisplayName', 'Overall DOS-wt (M1)');
    plot(E_avg_all_DOS_sca, sca_avg_all_DOS_tau_inv, 'b:', 'LineWidth', 3, ...
         'DisplayName', 'Overall DOS-wt (M2)');
end

% Format scattering rate plot
set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([0, 0.45]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('\tau^{-1} [fs^{-1}]', 'FontSize', 12);
title('ZrNiSn IIS Scattering Rates - EF=0, T=300K', 'FontSize', 14);

% Add legend for scattering rates
if ~isempty(legend_entries_sca)
    legend(legend_entries_sca, 'FontSize', 8, 'Location', 'best');
end

% Add annotation
text(0.02, 0.98, 'Method 1: DOS-wt S_k = \Sigma(S_k·DOS_k)/\Sigma(DOS_k)', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
text(0.02, 0.94, 'Method 2: 1/(\Sigma(\tau_k·DOS_k)/\Sigma(DOS_k))', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);

%% ===== FIGURE 2: RELAXATION TIMES =====
figure(figure2);

% Plot individual k-points as small dots for first band only
if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_tau_by_band{1}, 5, [0.7 0.7 0.7], 'filled', ...
            'DisplayName', 'Individual k-points (Band 1)');
end

% Plot per-band averages
h_tau = gobjects(numBands, 3);
legend_entries_tau = {};
legend_count = 0;

for b = 1:numBands
    if ~isempty(E_avg_by_band_simple_tau{b})
        % Simple average (solid line)
        h_tau(b,1) = plot(E_avg_by_band_simple_tau{b}, tau_avg_by_band_simple{b}, '-', ...
                         'Color', palette(b,:), 'LineWidth', 1.5);
        legend_count = legend_count + 1;
        legend_entries_tau{legend_count} = sprintf('Band %d - Simple', b);
        
        % DOS-weighted Method 1 (dashed line)
        h_tau(b,2) = plot(E_avg_by_band_DOS_tau{b}, tau_avg_by_band_DOS{b}, '--', ...
                         'Color', palette(b,:), 'LineWidth', 2);
        legend_count = legend_count + 1;
        legend_entries_tau{legend_count} = sprintf('Band %d - DOS-wt (Method 1)', b);
        
        % DOS-weighted Method 2 (dotted line)
        h_tau(b,3) = plot(E_avg_by_band_DOS_sca_inv{b}, tau_avg_by_band_DOS_sca_inv{b}, ':', ...
                         'Color', palette(b,:), 'LineWidth', 2);
        legend_count = legend_count + 1;
        legend_entries_tau{legend_count} = sprintf('Band %d - DOS-wt (Method 2)', b);
    end
end

% Plot overall averages
if ~isempty(E_avg_all_simple_tau)
    plot(E_avg_all_simple_tau, tau_avg_all_simple, 'k-', 'LineWidth', 3, ...
         'DisplayName', 'Overall Simple');
    plot(E_avg_all_DOS_tau, tau_avg_all_DOS, 'r--', 'LineWidth', 3, ...
         'DisplayName', 'Overall DOS-wt (M1)');
    plot(E_avg_all_DOS_tau, tau_avg_all_DOS_sca_inv, 'b:', 'LineWidth', 3, ...
         'DisplayName', 'Overall DOS-wt (M2)');
end

% Format relaxation time plot
set(gca, 'YScale', 'log');
ylim([1e-3, 1e1]);  % Adjust based on your data
xlim([0, 0.45]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('\tau [fs]', 'FontSize', 12);
title('ZrNiSn IIS Relaxation Times - EF=0, T=300K', 'FontSize', 14);

% Add legend for relaxation times
if ~isempty(legend_entries_tau)
    legend(legend_entries_tau, 'FontSize', 8, 'Location', 'best');
end

% Add annotation
text(0.02, 0.98, 'Method 1: DOS-wt \tau = \Sigma(\tau_k·DOS_k)/\Sigma(DOS_k)', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
text(0.02, 0.94, 'Method 2: 1/(\Sigma(S_k·DOS_k)/\Sigma(DOS_k))', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);

%% ===== SAVE DATA FOR XMGRACE =====
fprintf('\n=== Saving data for xmgrace ===\n');
baseName = 'ZrNiSn_IIS_300K';

% Create output directory
outputDir = 'IIS_xmgrace_300K_full';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% ===== SAVE SCATTERING RATE FILES =====
fprintf('\n--- Scattering Rate Files ---\n');

% Per-band scattering rates
for b = 1:numBands
    if ~isempty(E_avg_by_band_simple{b})
        % Simple average
        filename = fullfile(outputDir, sprintf('%s_band%02d_sca_simple.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - Simple Average Scattering Rate\n', b);
        fprintf(fid, '# Energy(eV)  Scattering_Rate(fs^-1)\n');
        for j = 1:length(E_avg_by_band_simple{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_simple{b}(j), sca_avg_by_band_simple{b}(j));
        end
        fclose(fid);
        
        % DOS-weighted Method 1
        filename = fullfile(outputDir, sprintf('%s_band%02d_sca_DOSwt_M1.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - DOS-weighted Scattering Rate (Method 1: weight S)\n', b);
        fprintf(fid, '# Energy(eV)  Scattering_Rate(fs^-1)\n');
        for j = 1:length(E_avg_by_band_DOS_sca{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_DOS_sca{b}(j), sca_avg_by_band_DOS{b}(j));
        end
        fclose(fid);
        
        % DOS-weighted Method 2
        filename = fullfile(outputDir, sprintf('%s_band%02d_sca_DOSwt_M2.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - DOS-weighted Scattering Rate (Method 2: weight tau, invert)\n', b);
        fprintf(fid, '# Energy(eV)  Scattering_Rate(fs^-1)\n');
        for j = 1:length(E_avg_by_band_DOS_tau_inv{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_DOS_tau_inv{b}(j), sca_avg_by_band_DOS_tau_inv{b}(j));
        end
        fclose(fid);
        
        fprintf('✓ Band %d scattering rates saved\n', b);
    end
end

% Overall scattering rates
filename = fullfile(outputDir, sprintf('%s_overall_sca_comparison.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# Overall Scattering Rate Comparison\n');
fprintf(fid, '# Energy  Simple  DOS-wt_M1  DOS-wt_M2  Diff_M1-M2  Ratio_M2/M1\n');
for j = 1:length(E_avg_all_simple_sca)
    diff = sca_avg_all_DOS_tau_inv(j) - sca_avg_all_DOS(j);
    ratio = sca_avg_all_DOS_tau_inv(j) / sca_avg_all_DOS(j);
    fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.4f\n', ...
            E_avg_all_simple_sca(j), sca_avg_all_simple(j), ...
            sca_avg_all_DOS(j), sca_avg_all_DOS_tau_inv(j), diff, ratio);
end
fclose(fid);

% ===== SAVE RELAXATION TIME FILES =====
fprintf('\n--- Relaxation Time Files ---\n');

% Per-band relaxation times
for b = 1:numBands
    if ~isempty(E_avg_by_band_simple_tau{b})
        % Simple average
        filename = fullfile(outputDir, sprintf('%s_band%02d_tau_simple.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - Simple Average Relaxation Time\n', b);
        fprintf(fid, '# Energy(eV)  Relaxation_Time(fs)\n');
        for j = 1:length(E_avg_by_band_simple_tau{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_simple_tau{b}(j), tau_avg_by_band_simple{b}(j));
        end
        fclose(fid);
        
        % DOS-weighted Method 1
        filename = fullfile(outputDir, sprintf('%s_band%02d_tau_DOSwt_M1.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - DOS-weighted Relaxation Time (Method 1: weight tau)\n', b);
        fprintf(fid, '# Energy(eV)  Relaxation_Time(fs)\n');
        for j = 1:length(E_avg_by_band_DOS_tau{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_DOS_tau{b}(j), tau_avg_by_band_DOS{b}(j));
        end
        fclose(fid);
        
        % DOS-weighted Method 2
        filename = fullfile(outputDir, sprintf('%s_band%02d_tau_DOSwt_M2.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - DOS-weighted Relaxation Time (Method 2: weight S, invert)\n', b);
        fprintf(fid, '# Energy(eV)  Relaxation_Time(fs)\n');
        for j = 1:length(E_avg_by_band_DOS_sca_inv{b})
            fprintf(fid, '%.6f  %.6e\n', E_avg_by_band_DOS_sca_inv{b}(j), tau_avg_by_band_DOS_sca_inv{b}(j));
        end
        fclose(fid);
        
        fprintf('✓ Band %d relaxation times saved\n', b);
    end
end

% Overall relaxation times
filename = fullfile(outputDir, sprintf('%s_overall_tau_comparison.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# Overall Relaxation Time Comparison\n');
fprintf(fid, '# Energy  Simple  DOS-wt_M1  DOS-wt_M2  Diff_M1-M2  Ratio_M2/M1\n');
for j = 1:length(E_avg_all_simple_tau)
    diff = tau_avg_all_DOS_sca_inv(j) - tau_avg_all_DOS(j);
    ratio = tau_avg_all_DOS_sca_inv(j) / tau_avg_all_DOS(j);
    fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.4f\n', ...
            E_avg_all_simple_tau(j), tau_avg_all_simple(j), ...
            tau_avg_all_DOS(j), tau_avg_all_DOS_sca_inv(j), diff, ratio);
end
fclose(fid);

% ===== SAVE COMBINED FILE WITH ALL DATA =====
fprintf('\n--- Combined file with all quantities ---\n');
filename = fullfile(outputDir, sprintf('%s_all_quantities.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# All quantities: scattering rates AND relaxation times\n');
fprintf(fid, '# Energy  S_simple  S_DOS_M1  S_DOS_M2  tau_simple  tau_DOS_M1  tau_DOS_M2\n');
for j = 1:min([length(E_avg_all_simple_sca), length(E_avg_all_DOS_sca), length(E_avg_all_DOS_tau)])
    fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
            E_avg_all_simple_sca(j), sca_avg_all_simple(j), ...
            sca_avg_all_DOS(j), sca_avg_all_DOS_tau_inv(j), ...
            tau_avg_all_simple(j), tau_avg_all_DOS(j), tau_avg_all_DOS_sca_inv(j));
end
fclose(fid);

% ===== SAVE INDIVIDUAL K-POINTS WITH BOTH SCATTERING AND RELAXATION =====
for b = 1:numBands
    if ~isempty(all_E_by_band{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_individual.dat', baseName, b));
        band_data = [all_E_by_band{b}, all_sca_by_band{b}, all_tau_by_band{b}, all_DOS_by_band{b}];
        dlmwrite(filename, band_data, 'delimiter', '\t', 'precision', '%.6e');
        fprintf('✓ Band %d individual data saved (E, S, tau, DOS)\n', b);
    end
end

% ===== CREATE XMGRACE BATCH FILES =====
fprintf('\n--- Creating xmgrace batch files ---\n');

% Batch file for scattering rates
filename = fullfile(outputDir, sprintf('%s_scattering_plot.bfile', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# xmgrace batch file for Scattering Rates\n');
fprintf(fid, 'title "ZrNiSn IIS Scattering Rates - EF=0, T=300K"\n');
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Scattering Rate [fs\\S-1\\N]"\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'legend on\n\n');
set_count = 0;
for b = 1:numBands
    if ~isempty(E_avg_by_band_simple{b})
        fprintf(fid, 'READ NXY "%s_band%02d_sca_simple.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d line style 1\n', set_count);
        fprintf(fid, 's%d legend "Band %d Simple"\n', set_count, b);
        set_count = set_count + 1;
        
        fprintf(fid, 'READ NXY "%s_band%02d_sca_DOSwt_M1.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d line style 2\n', set_count);
        fprintf(fid, 's%d legend "Band %d DOS-wt M1"\n', set_count, b);
        set_count = set_count + 1;
    end
end
fprintf(fid, 'READ NXY "%s_overall_sca_comparison.dat"\n', baseName);
fprintf(fid, 's%d line color 1\n', set_count);
fprintf(fid, 's%d line linewidth 3\n', set_count);
fprintf(fid, 'AUTOSCALE\n');
fclose(fid);

% Batch file for relaxation times
filename = fullfile(outputDir, sprintf('%s_relaxation_plot.bfile', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# xmgrace batch file for Relaxation Times\n');
fprintf(fid, 'title "ZrNiSn IIS Relaxation Times - EF=0, T=300K"\n');
fprintf(fid, 'xaxis label "Energy [eV]"\n');
fprintf(fid, 'yaxis label "Relaxation Time [fs]"\n');
fprintf(fid, 'yaxis scale logarithmic\n');
fprintf(fid, 'legend on\n\n');
set_count = 0;
for b = 1:numBands
    if ~isempty(E_avg_by_band_simple_tau{b})
        fprintf(fid, 'READ NXY "%s_band%02d_tau_simple.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d line style 1\n', set_count);
        fprintf(fid, 's%d legend "Band %d Simple"\n', set_count, b);
        set_count = set_count + 1;
        
        fprintf(fid, 'READ NXY "%s_band%02d_tau_DOSwt_M1.dat"\n', baseName, b);
        fprintf(fid, 's%d line color %d\n', set_count, b);
        fprintf(fid, 's%d line linewidth 2\n', set_count);
        fprintf(fid, 's%d line style 2\n', set_count);
        fprintf(fid, 's%d legend "Band %d DOS-wt M1"\n', set_count, b);
        set_count = set_count + 1;
    end
end
fprintf(fid, 'READ NXY "%s_overall_tau_comparison.dat"\n', baseName);
% Add these new arrays for the transport-weighted effective tau
% Initialize arrays for transport-weighted effective tau
E_eff_all_transport = [];
tau_eff_all_transport = [];

E_eff_by_band_transport = cell(numBands, 1);
tau_eff_by_band_transport = cell(numBands, 1);

% Configuration
EF_index = 8;            % Index for EF=0
T_index = 1;             % Temperature index (1=300K, adjust if needed)

fprintf('========================================\n');
fprintf('Processing IIS data for EF=0, T=300K\n');
fprintf('Calculating BOTH scattering rates AND relaxation times\n');
fprintf('Including transport-weighted effective tau\n');
fprintf('========================================\n');

%% Main processing loop
max_idx = 451;  % Fixed 451 energy points

for i = 1:max_idx
    % Arrays for this energy point (across bands)
    E_row = [];
    tau_eff_transport_row = [];
    DOS_row = [];  % For overall band weighting
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            % Get energy for this point
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= EF_index && size(tau_3D, 2) >= T_index
                % Extract tau(EF=0, T=300K, all k-points)
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                
                % Extract velocity squared (V^2) for all k-points
                % NOTE: You need to adjust this based on how velocities are stored
                if isfield(state_ID(i,b), 'velocity_squared') 
                    V2_k = state_ID(i,b).velocity_squared(:);
                elseif isfield(state_ID(i,b), 'vel')
                    % If velocity components are stored, calculate V^2
                    V2_k = state_ID(i,b).vel(:,1).^2 + state_ID(i,b).vel(:,2).^2 + state_ID(i,b).vel(:,3).^2;
                else
                    % If no velocity data, use group velocity from band structure
                    % This is a placeholder - adjust based on your data structure
                    warning('No velocity data at E=%.4f, band %d. Using uniform V^2.', E, b);
                    V2_k = ones(size(tau_k));
                end
                
                % Extract DOS for all k-points
                if isfield(state_ID(i,b), 'DOS') && ~isempty(state_ID(i,b).DOS)
                    DOS_k = state_ID(i,b).DOS(:);
                    
                    if length(DOS_k) ~= length(tau_k)
                        warning('DOS and tau size mismatch at E=%.4f, band %d.', E, b);
                        DOS_k = ones(size(tau_k));
                    end
                else
                    warning('No DOS data at E=%.4f, band %d. Using uniform weights.', E, b);
                    DOS_k = ones(size(tau_k));
                end
                
                % Filter: Keep only positive tau values
                positive_mask = tau_k > 0;
                
                if any(positive_mask)
                    tau_k_positive = tau_k(positive_mask);
                    V2_k_positive = V2_k(positive_mask);
                    DOS_k_positive = DOS_k(positive_mask);
                    sca_k = 1 ./ tau_k_positive;  % Scattering rates
                    
                    % ==== NEW: Calculate transport-weighted effective tau ====
                    % <tau> = <sum (1/S_k) * V^2_k * g_k>
                    % = sum( (1/S_k) * V^2_k * g_k ) / sum( V^2_k * g_k )
                    
                    % Weight = V^2 * g (transport weight)
                    transport_weight = V2_k_positive .* DOS_k_positive;
                    
                    % Effective tau using transport weight
                    if sum(transport_weight) > 0
                        tau_eff_band = sum(tau_k_positive .* transport_weight) / sum(transport_weight);
                    else
                        tau_eff_band = mean(tau_k_positive);  % Fallback to simple average
                    end
                    
                    % Store band-specific effective tau
                    E_eff_by_band_transport{b}(end+1) = E;
                    tau_eff_by_band_transport{b}(end+1) = tau_eff_band;
                    
                    % Collect for overall average
                    E_row(end+1) = E;
                    tau_eff_transport_row(end+1) = tau_eff_band;
                    DOS_row(end+1) = sum(DOS_k_positive);  % Total DOS for this band
                    
                    % ... [rest of your existing code remains the same] ...
                end
            end
        end
    end
    
    % Calculate overall transport-weighted effective tau
    if ~isempty(E_row) && ~isempty(DOS_row) && sum(DOS_row) > 0
        E_eff_all_transport(end+1) = sum(E_row .* DOS_row) / sum(DOS_row);
        tau_eff_all_transport(end+1) = sum(tau_eff_transport_row .* DOS_row) / sum(DOS_row);
    elseif ~isempty(E_row)
        E_eff_all_transport(end+1) = mean(E_row);
        tau_eff_all_transport(end+1) = mean(tau_eff_transport_row);
    end
    
    % Progress indicator
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', i, max_idx);
    end
end

%% ===== NEW FIGURE: TRANSPORT-WEIGHTED EFFECTIVE TAU =====
figure('Name', 'Transport-Weighted Effective Tau');
fig3 = gcf;
fig3.Position(3:4) = [1000, 800];

% Plot per-band effective tau
for b = 1:numBands
    if ~isempty(E_eff_by_band_transport{b})
        plot(E_eff_by_band_transport{b}, tau_eff_by_band_transport{b}, '-', ...
             'Color', palette(b,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('Band %d - Transport-wt', b));
        hold on;
    end
end

% Plot overall effective tau
if ~isempty(E_eff_all_transport)
    plot(E_eff_all_transport, tau_eff_all_transport, 'k-', 'LineWidth', 3, ...
         'DisplayName', 'Overall Transport-wt');
end

% Compare with other averaging methods
if ~isempty(E_avg_all_simple_tau)
    plot(E_avg_all_simple_tau, tau_avg_all_simple, 'r--', 'LineWidth', 1.5, ...
         'DisplayName', 'Overall Simple');
end
if ~isempty(E_avg_all_DOS_tau)
    plot(E_avg_all_DOS_tau, tau_avg_all_DOS, 'b--', 'LineWidth', 1.5, ...
         'DisplayName', 'Overall DOS-wt');
end

set(gca, 'YScale', 'log');
ylim([1e-3, 1e1]);
xlim([0, 0.45]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('\tau_{eff} [fs]', 'FontSize', 12);
title('Transport-Weighted Effective Relaxation Time', 'FontSize', 14);
legend('Location', 'best');

% Add annotation with formula
text(0.02, 0.98, '<\tau> = <\Sigma (1/S_k) V^2_k g_k>', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'tex');

%% ===== SAVE TRANSPORT-WEIGHTED DATA =====
fprintf('\n--- Saving transport-weighted effective tau ---\n');

% Per-band
for b = 1:numBands
    if ~isempty(E_eff_by_band_transport{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_tau_transport.dat', baseName, b));
        fid = fopen(filename, 'w');
        fprintf(fid, '# Band %d - Transport-weighted Effective Relaxation Time\n', b);
        fprintf(fid, '# Energy(eV)  Tau_effective(fs)\n');
        for j = 1:length(E_eff_by_band_transport{b})
            fprintf(fid, '%.6f  %.6e\n', E_eff_by_band_transport{b}(j), tau_eff_by_band_transport{b}(j));
        end
        fclose(fid);
        fprintf('✓ Band %d transport-weighted tau saved\n', b);
    end
end

% Overall comparison
filename = fullfile(outputDir, sprintf('%s_tau_methods_comparison.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# Comparison of Tau Averaging Methods\n');
fprintf(fid, '# Energy  Tau_simple  Tau_DOS-wt  Tau_transport-wt  Ratio_trans/DoS  Ratio_trans/Simple\n');
for j = 1:min([length(E_avg_all_simple_tau), length(E_avg_all_DOS_tau), length(E_eff_all_transport)])
    ratio_DOS = tau_eff_all_transport(j) / tau_avg_all_DOS(j);
    ratio_simple = tau_eff_all_transport(j) / tau_avg_all_simple(j);
    fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.4f  %.4f\n', ...
            E_avg_all_simple_tau(j), tau_avg_all_simple(j), ...
            tau_avg_all_DOS(j), tau_eff_all_transport(j), ...
            ratio_DOS, ratio_simple);
end
fclose(fid);

% Add statistics summary
fprintf('\n=== Transport-Weighted Tau Statistics ===\n');
fprintf('Mean ratio (Transport/DOS-wt): %.4f\n', mean(tau_eff_all_transport ./ tau_avg_all_DOS));
fprintf('Mean ratio (Transport/Simple): %.4f\n', mean(tau_eff_all_transport ./ tau_avg_all_simple));
fprintf('Max ratio (Transport/DOS-wt): %.4f\n', max(tau_eff_all_transport ./ tau_avg_all_DOS));
fprintf('Min ratio (Transport/Simple): %.4f\n', min(tau_eff_all_transport ./ tau_avg_all_simple));
