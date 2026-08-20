%% Load the data files
%load 'TE_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';

%% Initialize for plotting - Create figures
figure1 = figure('Name', 'Scattering Rates');
fig1 = gcf;
fig1.Position(3:4) = [1000, 800];
hold on;

figure2 = figure('Name', 'Relaxation Times');
fig2 = gcf;
fig2.Position(3:4) = [1000, 800];
hold on;

numBands = 3;
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

% Initialize arrays for data collection - FIVE METHODS
E_avg_all_M1 = []; sca_avg_all_M1 = []; tau_avg_all_M1 = [];
E_avg_all_M2 = []; sca_avg_all_M2 = []; tau_avg_all_M2 = [];
E_avg_all_M3 = []; sca_avg_all_M3 = []; tau_avg_all_M3 = [];
E_avg_all_M4 = []; sca_avg_all_M4 = []; tau_avg_all_M4 = [];
E_avg_all_M5 = []; sca_avg_all_M5 = []; tau_avg_all_M5 = [];

% Initialize per-band storage
E_avg_by_band_M1 = cell(numBands, 1); sca_avg_by_band_M1 = cell(numBands, 1); tau_avg_by_band_M1 = cell(numBands, 1);
E_avg_by_band_M2 = cell(numBands, 1); sca_avg_by_band_M2 = cell(numBands, 1); tau_avg_by_band_M2 = cell(numBands, 1);
E_avg_by_band_M3 = cell(numBands, 1); sca_avg_by_band_M3 = cell(numBands, 1); tau_avg_by_band_M3 = cell(numBands, 1);
E_avg_by_band_M4 = cell(numBands, 1); sca_avg_by_band_M4 = cell(numBands, 1); tau_avg_by_band_M4 = cell(numBands, 1);
E_avg_by_band_M5 = cell(numBands, 1); sca_avg_by_band_M5 = cell(numBands, 1); tau_avg_by_band_M5 = cell(numBands, 1);

% Store intermediate quantities for M5
E_avg_by_band_V2g = cell(numBands, 1);
V2g_avg_by_band = cell(numBands, 1);
E_avg_by_band_TDF = cell(numBands, 1);
TDF_avg_by_band = cell(numBands, 1);

% Individual k-points
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);
all_tau_by_band = cell(numBands, 1);
all_DOS_by_band = cell(numBands, 1);
all_V2_by_band = cell(numBands, 1);

% For combined storage
all_E = []; all_sca = []; all_tau = []; all_band = []; all_DOS = []; all_V2 = [];

% Configuration
EF_index = 8;            % Index for EF=0 (8th column in 2nd dimension)
T_index = 1;             % Temperature index (1st in 3rd dimension = 300K)

% Check for TDF variable
if exist('TDF', 'var')
    TDF_fields = fieldnames(TDF);
    fprintf('TDF struct found with fields: ');
    for f = 1:length(TDF_fields)
        fprintf('%s ', TDF_fields{f});
    end
    fprintf('\n');
    
    if isfield(TDF, 'xx')
        fprintf('Using TDF.xx for transport calculations\n');
        fprintf('TDF.xx dimensions: [%d, %d, %d]\n', size(TDF.xx));
    else
        error('TDF.xx field not found!');
    end
    use_TDF = true;
else
    fprintf('WARNING: TDF variable not found.\n');
    fprintf('M5 will use M4 as approximation.\n');
    use_TDF = false;
end

% Get energy range from E_array
if exist('E_array', 'var')
    E_min = min(E_array);
    E_max = max(E_array);
    fprintf('Using E_array: %.4f to %.4f eV (%d points)\n', E_min, E_max, length(E_array));
    
    energy_indices = [];
    for i = 1:size(state_ID, 1)
        if ~isempty(state_ID(i,1).E) && state_ID(i,1).E >= E_min && state_ID(i,1).E <= E_max
            energy_indices = [energy_indices, i];
        end
    end
    max_idx = length(energy_indices);
else
    max_idx = size(state_ID, 1);
    energy_indices = 1:max_idx;
    E_min = 0;
    E_max = 0.5;
    warning('E_array not found. Using all %d available energy points.', max_idx);
end

fprintf('========================================\n');
fprintf('Processing IIS data for EF=0, T=300K\n');
fprintf('Energy range: %.4f to %.4f eV\n', E_min, E_max);
fprintf('Number of energy points: %d\n', max_idx);
fprintf('EF_index = %d, T_index = %d\n', EF_index, T_index);
fprintf('========================================\n');
fprintf('\nFIVE AVERAGING METHODS:\n');
fprintf('M1: Simple <S>, tau = 1/<S>\n');
fprintf('M2: DOS-wt <S>, tau = 1/<S>\n');
fprintf('M3: DOS-wt <tau>, S = 1/<tau>\n');
fprintf('M4: Transport-wt <tau> (V^2*g), S = 1/<tau>\n');
fprintf('M5: TDF-based tau = TDF(E) / (V^2(E)*g(E)), S = 1/<tau>\n');

%% Main processing loop for M1-M4 and collecting V²·g for M5
for idx = 1:max_idx
    i = energy_indices(idx);
    
    E_row = [];
    sca_M1_row = []; tau_M1_row = [];
    sca_M2_row = []; tau_M2_row = [];
    sca_M3_row = []; tau_M3_row = [];
    sca_M4_row = []; tau_M4_row = [];
    DOS_row = [];
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= EF_index && size(tau_3D, 2) >= T_index
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                
                % Extract velocity
                if isfield(state_ID(i,b), 'V')
                    V_k = state_ID(i,b).V(:);
                    V2_k = V_k.^2;
                else
                    V2_k = ones(size(tau_k));
                end
                if length(V2_k) ~= length(tau_k), V2_k = ones(size(tau_k)); end
                
                % Extract DOS
                if isfield(state_ID(i,b), 'DOS') && ~isempty(state_ID(i,b).DOS)
                    DOS_k = state_ID(i,b).DOS(:);
                    if length(DOS_k) ~= length(tau_k), DOS_k = ones(size(tau_k)); end
                else
                    DOS_k = ones(size(tau_k));
                end
                
                positive_mask = tau_k > 0 & V2_k > 0;
                
                if any(positive_mask)
                    tau_k_positive = tau_k(positive_mask);
                    V2_k_positive = V2_k(positive_mask);
                    DOS_k_positive = DOS_k(positive_mask);
                    sca_k = 1 ./ tau_k_positive;
                    
                    % Store individual k-points
                    n_points = length(sca_k);
                    all_E_by_band{b} = [all_E_by_band{b}; repmat(E, n_points, 1)];
                    all_sca_by_band{b} = [all_sca_by_band{b}; sca_k];
                    all_tau_by_band{b} = [all_tau_by_band{b}; tau_k_positive];
                    all_DOS_by_band{b} = [all_DOS_by_band{b}; DOS_k_positive];
                    all_V2_by_band{b} = [all_V2_by_band{b}; V2_k_positive];
                    
                    all_E = [all_E; repmat(E, n_points, 1)];
                    all_sca = [all_sca; sca_k];
                    all_tau = [all_tau; tau_k_positive];
                    all_band = [all_band; repmat(b, n_points, 1)];
                    all_DOS = [all_DOS; DOS_k_positive];
                    all_V2 = [all_V2; V2_k_positive];
                    
                    % METHOD 1: Simple average of S
                    sca_M1 = mean(sca_k);
                    tau_M1 = 1 / sca_M1;
                    
                    % METHOD 2: DOS-weighted S
                    if sum(DOS_k_positive) > 0
                        sca_M2 = sum(sca_k .* DOS_k_positive) / sum(DOS_k_positive);
                    else
                        sca_M2 = mean(sca_k);
                    end
                    tau_M2 = 1 / sca_M2;
                    
                    % METHOD 3: DOS-weighted tau
                    if sum(DOS_k_positive) > 0
                        tau_M3 = sum(tau_k_positive .* DOS_k_positive) / sum(DOS_k_positive);
                    else
                        tau_M3 = mean(tau_k_positive);
                    end
                    sca_M3 = 1 / tau_M3;
                    
                    % METHOD 4: Transport-weighted tau (V²·g)
                    transport_weight_M4 = V2_k_positive .* DOS_k_positive;
                    if sum(transport_weight_M4) > 0
                        tau_M4 = sum(tau_k_positive .* transport_weight_M4) / sum(transport_weight_M4);
                    else
                        tau_M4 = mean(tau_k_positive);
                    end
                    sca_M4 = 1 / tau_M4;
                    
                    % Calculate V²·g sum for M5
                    V2g_sum = sum(V2_k_positive .* DOS_k_positive);
                    
                    % Store per-band results
                    E_avg_by_band_M1{b}(end+1) = E;
                    sca_avg_by_band_M1{b}(end+1) = sca_M1;
                    tau_avg_by_band_M1{b}(end+1) = tau_M1;
                    
                    E_avg_by_band_M2{b}(end+1) = E;
                    sca_avg_by_band_M2{b}(end+1) = sca_M2;
                    tau_avg_by_band_M2{b}(end+1) = tau_M2;
                    
                    E_avg_by_band_M3{b}(end+1) = E;
                    sca_avg_by_band_M3{b}(end+1) = sca_M3;
                    tau_avg_by_band_M3{b}(end+1) = tau_M3;
                    
                    E_avg_by_band_M4{b}(end+1) = E;
                    sca_avg_by_band_M4{b}(end+1) = sca_M4;
                    tau_avg_by_band_M4{b}(end+1) = tau_M4;
                    
                    % Store V²·g for M5
                    E_avg_by_band_V2g{b}(end+1) = E;
                    V2g_avg_by_band{b}(end+1) = V2g_sum;
                    
                    % Collect for overall averages
                    E_row(end+1) = E;
                    sca_M1_row(end+1) = sca_M1; tau_M1_row(end+1) = tau_M1;
                    sca_M2_row(end+1) = sca_M2; tau_M2_row(end+1) = tau_M2;
                    sca_M3_row(end+1) = sca_M3; tau_M3_row(end+1) = tau_M3;
                    sca_M4_row(end+1) = sca_M4; tau_M4_row(end+1) = tau_M4;
                    DOS_row(end+1) = sum(DOS_k_positive);
                end
            end
        end
    end
    
    % Overall averages across bands (M1-M4)
    if ~isempty(E_row) && sum(DOS_row) > 0
        E_avg = sum(E_row .* DOS_row) / sum(DOS_row);
        
        E_avg_all_M1(end+1) = E_avg;
        sca_avg_all_M1(end+1) = sum(sca_M1_row .* DOS_row) / sum(DOS_row);
        tau_avg_all_M1(end+1) = sum(tau_M1_row .* DOS_row) / sum(DOS_row);
        
        E_avg_all_M2(end+1) = E_avg;
        sca_avg_all_M2(end+1) = sum(sca_M2_row .* DOS_row) / sum(DOS_row);
        tau_avg_all_M2(end+1) = sum(tau_M2_row .* DOS_row) / sum(DOS_row);
        
        E_avg_all_M3(end+1) = E_avg;
        sca_avg_all_M3(end+1) = sum(sca_M3_row .* DOS_row) / sum(DOS_row);
        tau_avg_all_M3(end+1) = sum(tau_M3_row .* DOS_row) / sum(DOS_row);
        
        E_avg_all_M4(end+1) = E_avg;
        sca_avg_all_M4(end+1) = sum(sca_M4_row .* DOS_row) / sum(DOS_row);
        tau_avg_all_M4(end+1) = sum(tau_M4_row .* DOS_row) / sum(DOS_row);
    end
    
    if mod(idx, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', idx, max_idx);
    end
end

% Get actual number of overall energy points
N_overall = length(E_avg_all_M1);
fprintf('\nOverall energy points after processing: %d\n', N_overall);

%% METHOD 5: TDF-based tau = TDF.xx(E) / (V²(E)·g(E))
fprintf('\n=== Calculating M5: TDF-based tau = TDF.xx(E) / (V^2(E)*g(E)) ===\n');

if use_TDF && exist('TDF', 'var') && isfield(TDF, 'xx')
    % Extract TDF.xx at EF_index (8) and T_index (1)
    TDF_xx_extracted = TDF.xx(:, EF_index, T_index);
    
    % Initialize overall M5 arrays
    E_avg_all_M5 = zeros(N_overall, 1);
    tau_avg_all_M5 = zeros(N_overall, 1);
    sca_avg_all_M5 = zeros(N_overall, 1);
    
    % Per-band M5 calculation
    for b = 1:numBands
        n_band_points = length(E_avg_by_band_V2g{b});
        
        if n_band_points > 0
            E_avg_by_band_M5{b} = zeros(n_band_points, 1);
            tau_avg_by_band_M5{b} = zeros(n_band_points, 1);
            sca_avg_by_band_M5{b} = zeros(n_band_points, 1);
            
            for j = 1:n_band_points
                E_val = E_avg_by_band_V2g{b}(j);
                V2g_val = V2g_avg_by_band{b}(j);
                
                % Interpolate TDF.xx at this energy
                if exist('E_array', 'var') && ~isempty(E_array) && length(E_array) == size(TDF.xx, 1)
                    TDF_at_E = interp1(E_array, TDF_xx_extracted, E_val, 'linear', 0);
                else
                    % Use nearest index
                    [~, idx_e] = min(abs(E_array - E_val));
                    if idx_e <= length(TDF_xx_extracted)
                        TDF_at_E = TDF_xx_extracted(idx_e);
                    else
                        TDF_at_E = 0;
                    end
                end
                
                % Store TDF value
                E_avg_by_band_TDF{b}(j) = E_val;
                TDF_avg_by_band{b}(j) = TDF_at_E;
                
                % Calculate tau_M5 = TDF(E) / (V²(E)·g(E))
                if V2g_val > 0 && TDF_at_E > 0
                    tau_M5_val = TDF_at_E / V2g_val;
                else
                    tau_M5_val = NaN;
                end
                
                E_avg_by_band_M5{b}(j) = E_val;
                tau_avg_by_band_M5{b}(j) = tau_M5_val;
                if ~isnan(tau_M5_val) && tau_M5_val > 0
                    sca_avg_by_band_M5{b}(j) = 1 / tau_M5_val;
                else
                    sca_avg_by_band_M5{b}(j) = NaN;
                end
            end
        end
    end
    
    % Overall M5: Use the same energy points as M1-M4
    for j = 1:N_overall
        E_val = E_avg_all_M1(j);
        
        % Sum V²·g across bands at this energy point
        % We need to match the energy index between overall and per-band arrays
        V2g_total = 0;
        for b = 1:numBands
            % Find the closest energy in per-band V2g arrays
            if ~isempty(E_avg_by_band_V2g{b})
                [~, idx_b] = min(abs(E_avg_by_band_V2g{b} - E_val));
                if idx_b <= length(V2g_avg_by_band{b})
                    V2g_total = V2g_total + V2g_avg_by_band{b}(idx_b);
                end
            end
        end
        
        % Get TDF.xx at this energy
        if exist('E_array', 'var') && ~isempty(E_array) && length(E_array) == size(TDF.xx, 1)
            TDF_total = interp1(E_array, TDF_xx_extracted, E_val, 'linear', 0);
        else
            [~, idx_e] = min(abs(E_array - E_val));
            if idx_e <= length(TDF_xx_extracted)
                TDF_total = TDF_xx_extracted(idx_e);
            else
                TDF_total = 0;
            end
        end
        
        if V2g_total > 0 && TDF_total > 0
            tau_M5_overall = TDF_total / V2g_total;
        else
            tau_M5_overall = NaN;
        end
        
        E_avg_all_M5(j) = E_val;
        tau_avg_all_M5(j) = tau_M5_overall;
        if ~isnan(tau_M5_overall) && tau_M5_overall > 0
            sca_avg_all_M5(j) = 1 / tau_M5_overall;
        else
            sca_avg_all_M5(j) = NaN;
        end
    end
    
    % Remove NaN values from M5 arrays
    valid_M5 = ~isnan(tau_avg_all_M5) & tau_avg_all_M5 > 0;
    E_avg_all_M5_valid = E_avg_all_M5(valid_M5);
    tau_avg_all_M5_valid = tau_avg_all_M5(valid_M5);
    sca_avg_all_M5_valid = sca_avg_all_M5(valid_M5);
    
    fprintf('M5 calculation complete. %d valid points out of %d.\n', sum(valid_M5), N_overall);
else
    fprintf('WARNING: TDF.xx not available. M5 set equal to M4.\n');
    E_avg_all_M5_valid = E_avg_all_M4;
    tau_avg_all_M5_valid = tau_avg_all_M4;
    sca_avg_all_M5_valid = sca_avg_all_M4;
end

fprintf('\n=== Processing Complete ===\n');

%% FIGURE 1: SCATTERING RATES
figure(figure1);

if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_sca_by_band{1}, 5, [0.7 0.7 0.7], 'filled', ...
            'DisplayName', 'Individual k-points (Band 1)');
end

if ~isempty(E_avg_all_M1)
    plot(E_avg_all_M1, sca_avg_all_M1, 'k-', 'LineWidth', 2.5, 'DisplayName', 'M1: Simple <S>');
    plot(E_avg_all_M2, sca_avg_all_M2, 'b--', 'LineWidth', 2, 'DisplayName', 'M2: DOS-wt <S>');
    plot(E_avg_all_M3, sca_avg_all_M3, 'g-.', 'LineWidth', 2, 'DisplayName', 'M3: 1/<tau>_{DOS}');
    plot(E_avg_all_M4, sca_avg_all_M4, 'r:', 'LineWidth', 2.5, 'DisplayName', 'M4: 1/<tau>_{V2g}');
    if exist('E_avg_all_M5_valid', 'var')
        plot(E_avg_all_M5_valid, sca_avg_all_M5_valid, 'm-', 'LineWidth', 3, 'DisplayName', 'M5: 1/<tau>_{TDF}');
    end
end

set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([E_min, E_max]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Scattering Rate [fs^{-1}]', 'FontSize', 12);
title('ZrNiSn IIS - All Scattering Rate Methods', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 8);

%% FIGURE 2: RELAXATION TIMES
figure(figure2);

if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_tau_by_band{1}, 5, [0.7 0.7 0.7], 'filled', ...
            'DisplayName', 'Individual k-points (Band 1)');
end

if ~isempty(E_avg_all_M1)
    plot(E_avg_all_M1, tau_avg_all_M1, 'k-', 'LineWidth', 2.5, 'DisplayName', 'M1: 1/<S>');
    plot(E_avg_all_M2, tau_avg_all_M2, 'b--', 'LineWidth', 2, 'DisplayName', 'M2: 1/<S>_{DOS}');
    plot(E_avg_all_M3, tau_avg_all_M3, 'g-.', 'LineWidth', 2, 'DisplayName', 'M3: <tau>_{DOS}');
    plot(E_avg_all_M4, tau_avg_all_M4, 'r:', 'LineWidth', 2.5, 'DisplayName', 'M4: <tau>_{V2g}');
    if exist('E_avg_all_M5_valid', 'var')
        plot(E_avg_all_M5_valid, tau_avg_all_M5_valid, 'm-', 'LineWidth', 3, 'DisplayName', 'M5: TDF/(V^2g)');
    end
end

set(gca, 'YScale', 'log');
ylim([1e-3, 1e1]);
xlim([E_min, E_max]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Relaxation Time \tau [fs]', 'FontSize', 12);
title('ZrNiSn IIS - All Relaxation Time Methods', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 8);

%% STATISTICS
fprintf('\n========================================\n');
fprintf('COMPREHENSIVE METHOD COMPARISON\n');
fprintf('========================================\n');

methods_names = {'M1: Simple <S>', 'M2: DOS-wt <S>', 'M3: DOS-wt <tau>', ...
                 'M4: Transport (V^2g)', 'M5: TDF/(V^2g)'};

fprintf('\n--- Mean tau Values ---\n');
tau_M1_mean = mean(tau_avg_all_M1(~isnan(tau_avg_all_M1) & tau_avg_all_M1 > 0));
tau_M2_mean = mean(tau_avg_all_M2(~isnan(tau_avg_all_M2) & tau_avg_all_M2 > 0));
tau_M3_mean = mean(tau_avg_all_M3(~isnan(tau_avg_all_M3) & tau_avg_all_M3 > 0));
tau_M4_mean = mean(tau_avg_all_M4(~isnan(tau_avg_all_M4) & tau_avg_all_M4 > 0));

fprintf('M1: %.4f fs\n', tau_M1_mean);
fprintf('M2: %.4f fs\n', tau_M2_mean);
fprintf('M3: %.4f fs\n', tau_M3_mean);
fprintf('M4: %.4f fs\n', tau_M4_mean);

if exist('tau_avg_all_M5_valid', 'var')
    tau_M5_mean = mean(tau_avg_all_M5_valid);
    fprintf('M5: %.4f fs\n', tau_M5_mean);
    
    fprintf('\n--- Ratios Relative to M5 ---\n');
    % Interpolate M1-M4 to M5 energy grid for comparison
    tau_M1_interp = interp1(E_avg_all_M1, tau_avg_all_M1, E_avg_all_M5_valid, 'linear');
    tau_M2_interp = interp1(E_avg_all_M2, tau_avg_all_M2, E_avg_all_M5_valid, 'linear');
    tau_M3_interp = interp1(E_avg_all_M3, tau_avg_all_M3, E_avg_all_M5_valid, 'linear');
    tau_M4_interp = interp1(E_avg_all_M4, tau_avg_all_M4, E_avg_all_M5_valid, 'linear');
    
    fprintf('M1/M5: %.4f\n', nanmean(tau_M1_interp ./ tau_avg_all_M5_valid));
    fprintf('M2/M5: %.4f\n', nanmean(tau_M2_interp ./ tau_avg_all_M5_valid));
    fprintf('M3/M5: %.4f\n', nanmean(tau_M3_interp ./ tau_avg_all_M5_valid));
    fprintf('M4/M5: %.4f\n', nanmean(tau_M4_interp ./ tau_avg_all_M5_valid));
end

%% SAVE DATA
fprintf('\n=== Saving data ===\n');
baseName = 'ZrNiSn_IIS_300K_TDF_corrected';
outputDir = 'IIS_xmgrace_TDF_corrected';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

filename = fullfile(outputDir, sprintf('%s_overall_5methods.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# Complete 5-method comparison\n');
fprintf(fid, '# M5: tau = TDF.xx(E, EF=0, T=300K) / (V^2(E)*g(E))\n');
fprintf(fid, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5\n');

% Use M5 energy grid
if exist('E_avg_all_M5_valid', 'var')
    n_pts = length(E_avg_all_M5_valid);
    for j = 1:n_pts
        E_val = E_avg_all_M5_valid(j);
        tau1 = interp1(E_avg_all_M1, tau_avg_all_M1, E_val, 'linear', NaN);
        tau2 = interp1(E_avg_all_M2, tau_avg_all_M2, E_val, 'linear', NaN);
        tau3 = interp1(E_avg_all_M3, tau_avg_all_M3, E_val, 'linear', NaN);
        tau4 = interp1(E_avg_all_M4, tau_avg_all_M4, E_val, 'linear', NaN);
        tau5 = tau_avg_all_M5_valid(j);
        fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e\n', E_val, tau1, tau2, tau3, tau4, tau5);
    end
else
    n_pts = min([length(E_avg_all_M1), length(E_avg_all_M2), length(E_avg_all_M3), length(E_avg_all_M4)]);
    for j = 1:n_pts
        fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                E_avg_all_M1(j), tau_avg_all_M1(j), tau_avg_all_M2(j), ...
                tau_avg_all_M3(j), tau_avg_all_M4(j), tau_avg_all_M4(j));
    end
end
fclose(fid);

fprintf('✓ Data saved to: %s\n', outputDir);
fprintf('\n=== Complete ===\n');