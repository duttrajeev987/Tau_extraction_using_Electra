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
% ===== METHOD 1: Simple (average S, then invert) =====
E_avg_all_M1 = [];
sca_avg_all_M1 = [];
tau_avg_all_M1 = [];

% ===== METHOD 2: DOS-weighted S (average S with g, then invert) =====
E_avg_all_M2 = [];
sca_avg_all_M2 = [];
tau_avg_all_M2 = [];

% ===== METHOD 3: DOS-weighted τ (average τ directly with g) =====
E_avg_all_M3 = [];
sca_avg_all_M3 = [];
tau_avg_all_M3 = [];

% ===== METHOD 4: Transport-weighted τ (average τ with V²·g) =====
E_avg_all_M4 = [];
sca_avg_all_M4 = [];
tau_avg_all_M4 = [];

% ===== METHOD 5: TDF-weighted τ (average τ directly with TDF) =====
E_avg_all_M5 = [];
sca_avg_all_M5 = [];
tau_avg_all_M5 = [];

% Initialize per-band storage
E_avg_by_band_M1 = cell(numBands, 1);
sca_avg_by_band_M1 = cell(numBands, 1);
tau_avg_by_band_M1 = cell(numBands, 1);

E_avg_by_band_M2 = cell(numBands, 1);
sca_avg_by_band_M2 = cell(numBands, 1);
tau_avg_by_band_M2 = cell(numBands, 1);

E_avg_by_band_M3 = cell(numBands, 1);
sca_avg_by_band_M3 = cell(numBands, 1);
tau_avg_by_band_M3 = cell(numBands, 1);

E_avg_by_band_M4 = cell(numBands, 1);
sca_avg_by_band_M4 = cell(numBands, 1);
tau_avg_by_band_M4 = cell(numBands, 1);

E_avg_by_band_M5 = cell(numBands, 1);
sca_avg_by_band_M5 = cell(numBands, 1);
tau_avg_by_band_M5 = cell(numBands, 1);

% Store TDF values
E_avg_by_band_TDF = cell(numBands, 1);
TDF_avg_by_band = cell(numBands, 1);

% Individual k-points
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);
all_tau_by_band = cell(numBands, 1);
all_DOS_by_band = cell(numBands, 1);
all_V2_by_band = cell(numBands, 1);

% For combined storage
all_E = [];
all_sca = [];
all_tau = [];
all_band = [];
all_DOS = [];
all_V2 = [];
all_TDF_val = [];

% Configuration
EF_index = 8;            % Index for EF=0
T_index = 1;             % Temperature index (1=300K)

% Check for TDF variable
if exist('TDF', 'var')
    fprintf('TDF variable found. Dimensions: [%d, %d, %d]\n', size(TDF));
    fprintf('Will extract TDF at EF_index=%d, T_index=%d\n', EF_index, T_index);
    use_TDF = true;
else
    fprintf('WARNING: TDF variable not found in workspace.\n');
    fprintf('M5 will use V²·g as approximation for TDF.\n');
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
fprintf('M1: Simple <S>, τ = 1/<S>\n');
fprintf('M2: DOS-wt <S>, τ = 1/<S>\n');
fprintf('M3: DOS-wt <τ>, S = 1/<τ>\n');
fprintf('M4: Transport-wt <τ> (V²·g), S = 1/<τ>\n');
fprintf('M5: TDF-wt <τ> (TDF from file), S = 1/<τ>\n');

%% Main processing loop
for idx = 1:max_idx
    i = energy_indices(idx);
    
    E_row = [];
    sca_M1_row = []; tau_M1_row = [];
    sca_M2_row = []; tau_M2_row = [];
    sca_M3_row = []; tau_M3_row = [];
    sca_M4_row = []; tau_M4_row = [];
    sca_M5_row = []; tau_M5_row = [];
    DOS_row = [];
    TDF_band_row = [];
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= EF_index && size(tau_3D, 2) >= T_index
                % Extract tau(EF=0, T=300K, all k-points)
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                
                % Extract TDF for this (energy, EF, T)
                if use_TDF && exist('TDF', 'var')
                    if size(TDF, 1) >= i && size(TDF, 2) >= EF_index && size(TDF, 3) >= T_index
                        TDF_val = TDF(i, EF_index, T_index);
                        TDF_k = TDF_val * ones(size(tau_k));
                    else
                        if isfield(state_ID(i,b), 'V')
                            V_k = state_ID(i,b).V(:);
                            V2_k = V_k.^2;
                        else
                            V2_k = ones(size(tau_k));
                        end
                        if isfield(state_ID(i,b), 'DOS')
                            DOS_k = state_ID(i,b).DOS(:);
                        else
                            DOS_k = ones(size(tau_k));
                        end
                        TDF_k = V2_k .* DOS_k;
                        TDF_val = NaN;
                    end
                else
                    if isfield(state_ID(i,b), 'V')
                        V_k = state_ID(i,b).V(:);
                        V2_k = V_k.^2;
                    else
                        V2_k = ones(size(tau_k));
                    end
                    if isfield(state_ID(i,b), 'DOS')
                        DOS_k = state_ID(i,b).DOS(:);
                    else
                        DOS_k = ones(size(tau_k));
                    end
                    TDF_k = V2_k .* DOS_k;
                    TDF_val = NaN;
                end
                
                % Extract velocity
                if isfield(state_ID(i,b), 'V')
                    V_k = state_ID(i,b).V(:);
                    V2_k = V_k.^2;
                else
                    V2_k = ones(size(tau_k));
                end
                
                if length(V2_k) ~= length(tau_k)
                    V2_k = ones(size(tau_k));
                end
                
                % Extract DOS
                if isfield(state_ID(i,b), 'DOS') && ~isempty(state_ID(i,b).DOS)
                    DOS_k = state_ID(i,b).DOS(:);
                    if length(DOS_k) ~= length(tau_k)
                        DOS_k = ones(size(tau_k));
                    end
                else
                    DOS_k = ones(size(tau_k));
                end
                
                % Ensure TDF_k matches
                if length(TDF_k) ~= length(tau_k)
                    TDF_k = ones(size(tau_k));
                end
                
                % Filter positive values
                positive_mask = tau_k > 0 & V2_k > 0 & TDF_k >= 0;
                
                if any(positive_mask)
                    tau_k_positive = tau_k(positive_mask);
                    V2_k_positive = V2_k(positive_mask);
                    DOS_k_positive = DOS_k(positive_mask);
                    TDF_k_positive = TDF_k(positive_mask);
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
                    all_TDF_val = [all_TDF_val; repmat(TDF_val, n_points, 1)];
                    
                    % ===== METHOD 1: Simple average of S =====
                    sca_M1 = mean(sca_k);
                    tau_M1 = 1 / sca_M1;
                    
                    % ===== METHOD 2: DOS-weighted S =====
                    if sum(DOS_k_positive) > 0
                        sca_M2 = sum(sca_k .* DOS_k_positive) / sum(DOS_k_positive);
                    else
                        sca_M2 = mean(sca_k);
                    end
                    tau_M2 = 1 / sca_M2;
                    
                    % ===== METHOD 3: DOS-weighted τ =====
                    if sum(DOS_k_positive) > 0
                        tau_M3 = sum(tau_k_positive .* DOS_k_positive) / sum(DOS_k_positive);
                    else
                        tau_M3 = mean(tau_k_positive);
                    end
                    sca_M3 = 1 / tau_M3;
                    
                    % ===== METHOD 4: Transport-weighted τ (V²·g) =====
                    transport_weight_M4 = V2_k_positive .* DOS_k_positive;
                    if sum(transport_weight_M4) > 0
                        tau_M4 = sum(tau_k_positive .* transport_weight_M4) / sum(transport_weight_M4);
                    else
                        tau_M4 = mean(tau_k_positive);
                    end
                    sca_M4 = 1 / tau_M4;
                    
                    % ===== METHOD 5: TDF-weighted τ =====
                    if sum(TDF_k_positive) > 0
                        tau_M5 = sum(tau_k_positive .* TDF_k_positive) / sum(TDF_k_positive);
                    else
                        tau_M5 = mean(tau_k_positive);
                    end
                    sca_M5 = 1 / tau_M5;
                    
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
                    
                    E_avg_by_band_M5{b}(end+1) = E;
                    sca_avg_by_band_M5{b}(end+1) = sca_M5;
                    tau_avg_by_band_M5{b}(end+1) = tau_M5;
                    
                    % Store TDF values
                    E_avg_by_band_TDF{b}(end+1) = E;
                    TDF_avg_by_band{b}(end+1) = TDF_val;
                    
                    % Collect for overall averages
                    E_row(end+1) = E;
                    sca_M1_row(end+1) = sca_M1; tau_M1_row(end+1) = tau_M1;
                    sca_M2_row(end+1) = sca_M2; tau_M2_row(end+1) = tau_M2;
                    sca_M3_row(end+1) = sca_M3; tau_M3_row(end+1) = tau_M3;
                    sca_M4_row(end+1) = sca_M4; tau_M4_row(end+1) = tau_M4;
                    sca_M5_row(end+1) = sca_M5; tau_M5_row(end+1) = tau_M5;
                    DOS_row(end+1) = sum(DOS_k_positive);
                    TDF_band_row(end+1) = TDF_val;
                end
            end
        end
    end
    
    % Overall averages across bands
    if ~isempty(E_row)
        if sum(DOS_row) > 0
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
            
            % For M5, use TDF for band weighting
            if sum(TDF_band_row(~isnan(TDF_band_row))) > 0
                E_avg_all_M5(end+1) = sum(E_row(~isnan(TDF_band_row)) .* TDF_band_row(~isnan(TDF_band_row))) / sum(TDF_band_row(~isnan(TDF_band_row)));
                sca_avg_all_M5(end+1) = sum(sca_M5_row(~isnan(TDF_band_row)) .* TDF_band_row(~isnan(TDF_band_row))) / sum(TDF_band_row(~isnan(TDF_band_row)));
                tau_avg_all_M5(end+1) = sum(tau_M5_row(~isnan(TDF_band_row)) .* TDF_band_row(~isnan(TDF_band_row))) / sum(TDF_band_row(~isnan(TDF_band_row)));
            else
                E_avg_all_M5(end+1) = E_avg;
                sca_avg_all_M5(end+1) = sum(sca_M5_row .* DOS_row) / sum(DOS_row);
                tau_avg_all_M5(end+1) = sum(tau_M5_row .* DOS_row) / sum(DOS_row);
            end
        end
    end
    
    if mod(idx, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', idx, max_idx);
    end
end

fprintf('\n=== Processing Complete ===\n');

%% ===== FIGURE 1: SCATTERING RATES =====
figure(figure1);

if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_sca_by_band{1}, 5, [0.7 0.7 0.7], 'filled', ...
            'DisplayName', 'Individual k-points (Band 1)');
end

if ~isempty(E_avg_all_M1)
    plot(E_avg_all_M1, sca_avg_all_M1, 'k-', 'LineWidth', 2.5, 'DisplayName', 'M1: Simple <S>');
    plot(E_avg_all_M2, sca_avg_all_M2, 'b--', 'LineWidth', 2, 'DisplayName', 'M2: DOS-wt <S>');
    plot(E_avg_all_M3, sca_avg_all_M3, 'g-.', 'LineWidth', 2, 'DisplayName', 'M3: 1/<τ>_{DOS}');
    plot(E_avg_all_M4, sca_avg_all_M4, 'r:', 'LineWidth', 2.5, 'DisplayName', 'M4: 1/<τ>_{V²g}');
    plot(E_avg_all_M5, sca_avg_all_M5, 'm-', 'LineWidth', 3, 'DisplayName', 'M5: 1/<τ>_{TDF}');
end

set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([E_min, E_max]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Scattering Rate [fs^{-1}]', 'FontSize', 12);
title('ZrNiSn IIS - All Scattering Rate Methods', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 8);

%% ===== FIGURE 2: RELAXATION TIMES =====
figure(figure2);

if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_tau_by_band{1}, 5, [0.7 0.7 0.7], 'filled', ...
            'DisplayName', 'Individual k-points (Band 1)');
end

if ~isempty(E_avg_all_M1)
    plot(E_avg_all_M1, tau_avg_all_M1, 'k-', 'LineWidth', 2.5, 'DisplayName', 'M1: 1/<S>');
    plot(E_avg_all_M2, tau_avg_all_M2, 'b--', 'LineWidth', 2, 'DisplayName', 'M2: 1/<S>_{DOS}');
    plot(E_avg_all_M3, tau_avg_all_M3, 'g-.', 'LineWidth', 2, 'DisplayName', 'M3: <τ>_{DOS}');
    plot(E_avg_all_M4, tau_avg_all_M4, 'r:', 'LineWidth', 2.5, 'DisplayName', 'M4: <τ>_{V²g}');
    plot(E_avg_all_M5, tau_avg_all_M5, 'm-', 'LineWidth', 3, 'DisplayName', 'M5: <τ>_{TDF}');
end

set(gca, 'YScale', 'log');
ylim([1e-3, 1e1]);
xlim([E_min, E_max]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Relaxation Time \tau [fs]', 'FontSize', 12);
title('ZrNiSn IIS - All Relaxation Time Methods', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 8);

%% ===== FIGURE 3: COMPREHENSIVE METHOD COMPARISON =====
figure('Name', 'Method Comparisons', 'Position', [100, 100, 1400, 1000]);

% Subplot 1: All taus together
subplot(3,4,1);
if ~isempty(E_avg_all_M1)
    plot(E_avg_all_M1, tau_avg_all_M1, 'k-', 'LineWidth', 2); hold on;
    plot(E_avg_all_M2, tau_avg_all_M2, 'b--', 'LineWidth', 2);
    plot(E_avg_all_M3, tau_avg_all_M3, 'g-.', 'LineWidth', 2);
    plot(E_avg_all_M4, tau_avg_all_M4, 'r:', 'LineWidth', 2);
    plot(E_avg_all_M5, tau_avg_all_M5, 'm-', 'LineWidth', 2);
end
set(gca, 'YScale', 'log');
xlabel('Energy [eV]'); ylabel('\tau [fs]');
title('All Methods - \tau');
legend('M1', 'M2', 'M3', 'M4', 'M5', 'Location', 'best', 'FontSize', 7);
grid on;

% Subplot 2: M1 vs M2 (Harmonic means)
subplot(3,4,2);
if ~isempty(tau_avg_all_M1)
    plot(E_avg_all_M1, tau_avg_all_M1, 'k-', 'LineWidth', 2); hold on;
    plot(E_avg_all_M2, tau_avg_all_M2, 'b--', 'LineWidth', 2);
    ratio12 = tau_avg_all_M2 ./ tau_avg_all_M1;
end
set(gca, 'YScale', 'log');
xlabel('Energy [eV]'); ylabel('\tau [fs]');
title(sprintf('M1 vs M2 (Harmonic)\nM2/M1 = %.4f', mean(ratio12)));
legend('M1: Simple', 'M2: DOS-wt S', 'Location', 'best', 'FontSize', 7);
grid on;

% Subplot 3: M3 vs M4 vs M5 (Arithmetic means)
subplot(3,4,3);
if ~isempty(tau_avg_all_M3)
    plot(E_avg_all_M3, tau_avg_all_M3, 'g-.', 'LineWidth', 2); hold on;
    plot(E_avg_all_M4, tau_avg_all_M4, 'r:', 'LineWidth', 2);
    plot(E_avg_all_M5, tau_avg_all_M5, 'm-', 'LineWidth', 2);
end
set(gca, 'YScale', 'log');
xlabel('Energy [eV]'); ylabel('\tau [fs]');
title('M3 vs M4 vs M5 (Arithmetic)');
legend('M3: DOS', 'M4: V²g', 'M5: TDF', 'Location', 'best', 'FontSize', 7);
grid on;

% Subplot 4: Harmonic vs Arithmetic (M2 vs M3)
subplot(3,4,4);
if ~isempty(tau_avg_all_M2)
    plot(E_avg_all_M2, tau_avg_all_M2, 'b--', 'LineWidth', 2); hold on;
    plot(E_avg_all_M3, tau_avg_all_M3, 'g-.', 'LineWidth', 2);
    ratio23 = tau_avg_all_M3 ./ tau_avg_all_M2;
end
set(gca, 'YScale', 'log');
xlabel('Energy [eV]'); ylabel('\tau [fs]');
title(sprintf('Harmonic vs Arithmetic\nM3/M2 = %.4f', mean(ratio23)));
legend('M2: Harm', 'M3: Arith', 'Location', 'best', 'FontSize', 7);
grid on;

% Subplot 5: All ratios vs M5
subplot(3,4,5);
if ~isempty(E_avg_all_M5)
    plot(E_avg_all_M5, tau_avg_all_M1 ./ tau_avg_all_M5, 'k-', 'LineWidth', 1.5); hold on;
    plot(E_avg_all_M5, tau_avg_all_M2 ./ tau_avg_all_M5, 'b--', 'LineWidth', 1.5);
    plot(E_avg_all_M5, tau_avg_all_M3 ./ tau_avg_all_M5, 'g-.', 'LineWidth', 1.5);
    plot(E_avg_all_M5, tau_avg_all_M4 ./ tau_avg_all_M5, 'r:', 'LineWidth', 1.5);
end
xlabel('Energy [eV]'); ylabel('Ratio to M5');
title('All Methods / M5 (TDF)');
legend('M1/M5', 'M2/M5', 'M3/M5', 'M4/M5', 'Location', 'best', 'FontSize', 7);
grid on;
yline(1, 'k--');

% Subplot 6: M4/M5 ratio
subplot(3,4,6);
if ~isempty(E_avg_all_M5)
    ratio45 = tau_avg_all_M4 ./ tau_avg_all_M5;
    plot(E_avg_all_M5, ratio45, 'm-', 'LineWidth', 2);
end
xlabel('Energy [eV]'); ylabel('M4/M5 Ratio');
title(sprintf('V²·g / TDF\nMean = %.4f, Std = %.4f', mean(ratio45), std(ratio45)));
grid on;
yline(1, 'k--');

% Subplot 7: M3/M1 ratio (Arith/Harm with DOS)
subplot(3,4,7);
if ~isempty(E_avg_all_M3)
    ratio31 = tau_avg_all_M3 ./ tau_avg_all_M1;
    plot(E_avg_all_M3, ratio31, 'c-', 'LineWidth', 2);
end
xlabel('Energy [eV]'); ylabel('M3/M1 Ratio');
title(sprintf('DOS-Arith / Simple-Harm\nMean = %.4f', mean(ratio31)));
grid on;
yline(1, 'k--');

% Subplot 8: Scattering rates M1 vs M5
subplot(3,4,8);
if ~isempty(E_avg_all_M1)
    loglog(E_avg_all_M1, sca_avg_all_M1, 'k-', 'LineWidth', 2); hold on;
    loglog(E_avg_all_M5, sca_avg_all_M5, 'm-', 'LineWidth', 2);
end
xlabel('Energy [eV]'); ylabel('S [fs^{-1}]');
title('Scattering: M1 vs M5');
legend('M1: Simple', 'M5: TDF', 'Location', 'best', 'FontSize', 7);
grid on;

% Subplot 9-11: Per-band velocity-scattering correlation
for b = 1:numBands
    subplot(3,4,8+b);
    if ~isempty(all_V2_by_band{b}) && ~isempty(all_sca_by_band{b})
        loglog(all_V2_by_band{b}, all_sca_by_band{b}, '.', 'Color', palette(b,:), 'MarkerSize', 6);
        xlabel('V^2 [(m/s)^2]', 'FontSize', 9);
        ylabel('S [fs^{-1}]', 'FontSize', 9);
        corr_val = corr(all_V2_by_band{b}, all_sca_by_band{b});
        
        % Power law fit
        if length(all_V2_by_band{b}) > 10
            valid = all_V2_by_band{b} > 0 & all_sca_by_band{b} > 0;
            if sum(valid) > 10
                p = polyfit(log10(all_V2_by_band{b}(valid)), log10(all_sca_by_band{b}(valid)), 1);
                hold on;
                x_fit = logspace(log10(min(all_V2_by_band{b}(valid))), log10(max(all_V2_by_band{b}(valid))), 100);
                y_fit = 10.^polyval(p, log10(x_fit));
                loglog(x_fit, y_fit, 'k-', 'LineWidth', 1.5);
            end
        end
        title(sprintf('Band %d: Corr=%.3f', b, corr_val), 'FontSize', 9);
        grid on;
    end
end

% Subplot 12: Summary statistics
subplot(3,4,12);
axis off;
text_str = sprintf('MEAN τ [fs]:\n');
tau_means = [mean(tau_avg_all_M1), mean(tau_avg_all_M2), mean(tau_avg_all_M3), ...
             mean(tau_avg_all_M4), mean(tau_avg_all_M5)];
for m = 1:5
    text_str = [text_str sprintf('M%d: %.4f\n', m, tau_means(m))];
end
text_str = [text_str sprintf('\nRATIOS vs M5:\n')];
for m = 1:4
    text_str = [text_str sprintf('M%d/M5: %.4f\n', m, tau_means(m)/tau_means(5))];
end
text(0.1, 0.5, text_str, 'FontSize', 8, 'VerticalAlignment', 'middle');

sgtitle('Comprehensive Method Comparison - ZrNiSn IIS', 'FontSize', 14, 'FontWeight', 'bold');

%% ===== FIGURE 4: VELOCITY-SCATTERING DETAILED ANALYSIS =====
figure('Name', 'Velocity-Scattering Analysis', 'Position', [150, 150, 1200, 800]);

for b = 1:numBands
    % Velocity-Scattering scatter with DOS color
    subplot(2, 3, b);
    if ~isempty(all_V2_by_band{b}) && ~isempty(all_sca_by_band{b}) && ~isempty(all_DOS_by_band{b})
        scatter(log10(all_V2_by_band{b}), log10(all_sca_by_band{b}), 8, ...
               all_DOS_by_band{b}, 'filled');
        colorbar;
        xlabel('log_{10}(V^2)'); ylabel('log_{10}(S)');
        corr_val = corr(all_V2_by_band{b}, all_sca_by_band{b});
        title(sprintf('Band %d: Corr(V²,S)=%.4f', b, corr_val));
        grid on;
        
        % Power law fit
        if length(all_V2_by_band{b}) > 10
            valid = all_V2_by_band{b} > 0 & all_sca_by_band{b} > 0;
            if sum(valid) > 10
                p = polyfit(log10(all_V2_by_band{b}(valid)), log10(all_sca_by_band{b}(valid)), 1);
                hold on;
                x_fit = linspace(min(log10(all_V2_by_band{b}(valid))), max(log10(all_V2_by_band{b}(valid))), 100);
                y_fit = polyval(p, x_fit);
                plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
                text(0.05, 0.95, sprintf('S ∝ (V^2)^{%.3f}', p(1)), ...
                     'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
            end
        end
    end
    
    % Velocity histogram
    subplot(2, 3, b+3);
    if ~isempty(all_V2_by_band{b})
        histogram(log10(all_V2_by_band{b}), 50, 'FaceColor', palette(b,:), 'EdgeColor', 'none');
        xlabel('log_{10}(V^2)'); ylabel('Count');
        title(sprintf('Band %d: V² Distribution\nCV=%.1f%%', b, 100*std(all_V2_by_band{b})/mean(all_V2_by_band{b})));
        grid on;
    end
end
sgtitle('Velocity-Scattering Rate Correlation Analysis', 'FontSize', 14, 'FontWeight', 'bold');

%% ===== FIGURE 5: TDF ANALYSIS =====
figure('Name', 'TDF Analysis', 'Position', [200, 200, 1200, 600]);

% TDF vs Energy
subplot(2,3,1);
for b = 1:numBands
    if ~isempty(E_avg_by_band_TDF{b}) && ~isempty(TDF_avg_by_band{b})
        semilogy(E_avg_by_band_TDF{b}, TDF_avg_by_band{b}, '-', ...
                 'Color', palette(b,:), 'LineWidth', 2);
        hold on;
    end
end
xlabel('Energy [eV]'); ylabel('TDF');
title('TDF(E, EF=0, T=300K)');
legend('Band 1', 'Band 2', 'Band 3', 'Location', 'best');
grid on;
xlim([E_min, E_max]);

% TDF ratio between bands
subplot(2,3,2);
if ~isempty(TDF_avg_by_band{1}) && ~isempty(TDF_avg_by_band{2})
    ratio_TDF_12 = TDF_avg_by_band{1} ./ TDF_avg_by_band{2};
    semilogy(E_avg_by_band_TDF{1}, ratio_TDF_12, 'b-', 'LineWidth', 2); hold on;
end
if ~isempty(TDF_avg_by_band{1}) && ~isempty(TDF_avg_by_band{3})
    ratio_TDF_13 = TDF_avg_by_band{1} ./ TDF_avg_by_band{3};
    semilogy(E_avg_by_band_TDF{1}, ratio_TDF_13, 'r-', 'LineWidth', 2);
end
xlabel('Energy [eV]'); ylabel('TDF Ratio');
title('TDF Band Ratios');
legend('Band1/Band2', 'Band1/Band3', 'Location', 'best');
grid on;
yline(1, 'k--');

% M5 vs M4 comparison
subplot(2,3,3);
if ~isempty(E_avg_all_M4)
    plot(E_avg_all_M4, tau_avg_all_M4, 'r:', 'LineWidth', 2.5); hold on;
    plot(E_avg_all_M5, tau_avg_all_M5, 'm-', 'LineWidth', 2.5);
end
set(gca, 'YScale', 'log');
xlabel('Energy [eV]'); ylabel('\tau [fs]');
title('M4 (V²·g) vs M5 (TDF)');
legend('M4', 'M5', 'Location', 'best');
grid on;

% M4/M5 ratio vs Energy
subplot(2,3,4);
if ~isempty(E_avg_all_M5)
    ratio45 = tau_avg_all_M4 ./ tau_avg_all_M5;
    plot(E_avg_all_M5, ratio45, 'm-', 'LineWidth', 2);
end
xlabel('Energy [eV]'); ylabel('M4/M5 Ratio');
title(sprintf('V²·g / TDF Ratio\nMean=%.4f, Std=%.4f', mean(ratio45), std(ratio45)));
grid on;
yline(1, 'k--');

% TDF vs V²·g correlation
subplot(2,3,5);
if ~isempty(all_V2) && ~isempty(all_DOS) && ~isempty(all_TDF_val)
    V2g_all = all_V2 .* all_DOS;
    % Remove NaN values
    valid_idx = ~isnan(all_TDF_val);
    if sum(valid_idx) > 0
        loglog(V2g_all(valid_idx), all_TDF_val(valid_idx), '.', 'Color', palette(1,:));
        xlabel('V²·g'); ylabel('TDF');
        corr_tdf = corr(V2g_all(valid_idx), all_TDF_val(valid_idx));
        title(sprintf('TDF vs V²·g\nCorr=%.4f', corr_tdf));
        grid on;
    end
end

% Summary text
subplot(2,3,6);
axis off;
text_str = sprintf('TDF ANALYSIS SUMMARY:\n\n');
if exist('ratio45', 'var')
    text_str = [text_str sprintf('M4/M5 mean ratio: %.4f\n', mean(ratio45))];
    if abs(mean(ratio45)-1) < 0.02
        text_str = [text_str sprintf('→ V²·g is EXCELLENT approx.\n')];
    elseif abs(mean(ratio45)-1) < 0.05
        text_str = [text_str sprintf('→ V²·g is GOOD approx.\n')];
    else
        text_str = [text_str sprintf('→ V²·g differs significantly\n')];
    end
end
text_str = [text_str sprintf('\nTDF includes full transport\nweighting (Fermi window etc.)\n')];
text_str = [text_str sprintf('\nM5 (TDF) is MOST ACCURATE\nfor transport calculations')];
text(0.1, 0.5, text_str, 'FontSize', 10, 'VerticalAlignment', 'middle');

sgtitle('TDF Analysis - Transport Distribution Function', 'FontSize', 14, 'FontWeight', 'bold');

%% ===== COMPREHENSIVE STATISTICS =====
fprintf('\n========================================\n');
fprintf('COMPREHENSIVE METHOD COMPARISON\n');
fprintf('========================================\n');

methods_names = {'M1: Simple <S>', 'M2: DOS-wt <S>', 'M3: DOS-wt <τ>', ...
                 'M4: Transport (V²·g)', 'M5: TDF-weighted'};
tau_arrays = {tau_avg_all_M1, tau_avg_all_M2, tau_avg_all_M3, ...
              tau_avg_all_M4, tau_avg_all_M5};
sca_arrays = {sca_avg_all_M1, sca_avg_all_M2, sca_avg_all_M3, ...
              sca_avg_all_M4, sca_avg_all_M5};

fprintf('\n--- Mean τ Values ---\n');
for m = 1:5
    fprintf('%-25s: %.4f fs\n', methods_names{m}, mean(tau_arrays{m}));
end

fprintf('\n--- Mean Scattering Rates ---\n');
for m = 1:5
    fprintf('%-25s: %.4e fs^-1\n', methods_names{m}, mean(sca_arrays{m}));
end

fprintf('\n--- Ratios Relative to M5 (TDF - Most Accurate) ---\n');
for m = 1:4
    ratio = tau_arrays{m} ./ tau_arrays{5};
    fprintf('%-25s: mean=%.4f, min=%.4f, max=%.4f, std=%.4f\n', ...
            [methods_names{m} '/M5'], mean(ratio), min(ratio), max(ratio), std(ratio));
end

fprintf('\n--- M4 vs M5 Comparison (V²·g vs TDF) ---\n');
if ~isempty(tau_avg_all_M4) && ~isempty(tau_avg_all_M5)
    ratio45 = tau_avg_all_M4 ./ tau_avg_all_M5;
    fprintf('M4/M5 ratio: mean=%.4f, min=%.4f, max=%.4f, std=%.4f\n', ...
            mean(ratio45), min(ratio45), max(ratio45), std(ratio45));
    fprintf('Percentage difference: %.2f%%\n', 100*abs(mean(ratio45)-1));
end

fprintf('\n--- Harmonic vs Arithmetic (M2 vs M3) ---\n');
if ~isempty(tau_avg_all_M2) && ~isempty(tau_avg_all_M3)
    ratio23 = tau_avg_all_M2 ./ tau_avg_all_M3;
    fprintf('M2/M3 (Harmonic/Arithmetic): mean=%.4f\n', mean(ratio23));
end

%% ===== VELOCITY ANALYSIS =====
fprintf('\n========================================\n');
fprintf('VELOCITY-SCATTERING ANALYSIS\n');
fprintf('========================================\n');

for b = 1:numBands
    if ~isempty(all_V2_by_band{b})
        V2_band = all_V2_by_band{b};
        fprintf('\nBand %d:\n', b);
        fprintf('  V² range: %.2e to %.2e\n', min(V2_band), max(V2_band));
        fprintf('  V² mean: %.2e ± %.2e\n', mean(V2_band), std(V2_band));
        fprintf('  V² CV: %.2f%%\n', 100*std(V2_band)/mean(V2_band));
        
        if ~isempty(all_sca_by_band{b})
            sca_band = all_sca_by_band{b};
            corr_val = corr(V2_band, sca_band);
            fprintf('  Correlation(V², S): %.4f\n', corr_val);
            
            % Power law fit
            valid = V2_band > 0 & sca_band > 0;
            if sum(valid) > 10
                p = polyfit(log10(V2_band(valid)), log10(sca_band(valid)), 1);
                fprintf('  Power law: S ∝ (V²)^{%.3f}\n', p(1));
                
                if p(1) < -0.3
                    fprintf('  → Strong IIS character (faster = less scattering)\n');
                elseif p(1) < -0.1
                    fprintf('  → Moderate IIS character\n');
                elseif abs(p(1)) < 0.1
                    fprintf('  → Velocity-independent scattering\n');
                else
                    fprintf('  → Non-IIS: faster = more scattering\n');
                end
            end
        end
    end
end

%% ===== TDF ANALYSIS =====
fprintf('\n========================================\n');
fprintf('TDF ANALYSIS\n');
fprintf('========================================\n');

if ~isempty(all_TDF_val) && ~isempty(all_V2) && ~isempty(all_DOS)
    V2g_all = all_V2 .* all_DOS;
    valid_idx = ~isnan(all_TDF_val);
    if sum(valid_idx) > 0
        corr_tdf = corr(V2g_all(valid_idx), all_TDF_val(valid_idx));
        fprintf('Correlation(TDF, V²·g): %.4f\n', corr_tdf);
        if corr_tdf > 0.95
            fprintf('→ V²·g is EXCELLENT approximation for TDF\n');
        elseif corr_tdf > 0.8
            fprintf('→ V²·g is GOOD approximation for TDF\n');
        else
            fprintf('→ V²·g differs from TDF significantly\n');
            fprintf('  Use TDF (M5) for accurate transport calculations\n');
        end
    end
end

%% ===== PHYSICAL RECOMMENDATIONS =====
fprintf('\n========================================\n');
fprintf('PHYSICAL RECOMMENDATIONS\n');
fprintf('========================================\n');

fprintf('\n1. BEST METHOD FOR TRANSPORT: M5 (TDF-weighted)\n');
fprintf('   τ_TDF = %.4f fs\n', mean(tau_avg_all_M5));

fprintf('\n2. If TDF unavailable: M4 (V²·g-weighted)\n');
fprintf('   τ_V²g = %.4f fs\n', mean(tau_avg_all_M4));

fprintf('\n3. For scattering analysis: M2 (DOS-wt harmonic)\n');
fprintf('   τ_harm = %.4f fs\n', mean(tau_avg_all_M2));

fprintf('\n4. Error using simpler methods vs M5:\n');
for m = 1:3
    error_pct = 100 * (mean(tau_arrays{m}) / mean(tau_avg_all_M5) - 1);
    fprintf('   %s: %+.1f%% error\n', methods_names{m}, error_pct);
end

%% ===== SAVE DATA =====
fprintf('\n=== Saving data ===\n');
baseName = 'ZrNiSn_IIS_300K_5methods_TDF';
outputDir = 'IIS_xmgrace_5methods_TDF_complete';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save overall comparison
filename = fullfile(outputDir, sprintf('%s_overall_5methods.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# Complete 5-method comparison with TDF\n');
fprintf(fid, '# TDF extracted at EF_index=%d, T_index=%d\n', EF_index, T_index);
fprintf(fid, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5');
fprintf(fid, '  S_M1  S_M2  S_M3  S_M4  S_M5');
fprintf(fid, '  Ratio_M1M5  Ratio_M2M5  Ratio_M3M5  Ratio_M4M5\n');
for j = 1:min([length(E_avg_all_M1), length(E_avg_all_M5)])
    fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.4f  %.4f  %.4f  %.4f\n', ...
            E_avg_all_M1(j), ...
            tau_avg_all_M1(j), tau_avg_all_M2(j), tau_avg_all_M3(j), ...
            tau_avg_all_M4(j), tau_avg_all_M5(j), ...
            sca_avg_all_M1(j), sca_avg_all_M2(j), sca_avg_all_M3(j), ...
            sca_avg_all_M4(j), sca_avg_all_M5(j), ...
            tau_avg_all_M1(j)/tau_avg_all_M5(j), tau_avg_all_M2(j)/tau_avg_all_M5(j), ...
            tau_avg_all_M3(j)/tau_avg_all_M5(j), tau_avg_all_M4(j)/tau_avg_all_M5(j));
end
fclose(fid);

% Save TDF values
if ~isempty(E_avg_by_band_TDF{1})
    filename = fullfile(outputDir, sprintf('%s_TDF_values.dat', baseName));
    fid = fopen(filename, 'w');
    fprintf(fid, '# TDF values at EF=0, T=300K\n');
    fprintf(fid, '# Energy  TDF_Band1  TDF_Band2  TDF_Band3\n');
    max_len = max([length(E_avg_by_band_TDF{1}), length(E_avg_by_band_TDF{2}), length(E_avg_by_band_TDF{3})]);
    for j = 1:max_len
        e = NaN; t1 = NaN; t2 = NaN; t3 = NaN;
        if j <= length(E_avg_by_band_TDF{1})
            e = E_avg_by_band_TDF{1}(j);
            t1 = TDF_avg_by_band{1}(j);
        end
        if j <= length(E_avg_by_band_TDF{2})
            t2 = TDF_avg_by_band{2}(j);
        end
        if j <= length(E_avg_by_band_TDF{3})
            t3 = TDF_avg_by_band{3}(j);
        end
        fprintf(fid, '%.6f  %.6e  %.6e  %.6e\n', e, t1, t2, t3);
    end
    fclose(fid);
end

% Save per-band individual k-points
for b = 1:numBands
    if ~isempty(all_E_by_band{b})
        filename = fullfile(outputDir, sprintf('%s_band%02d_kpoints.dat', baseName, b));
        data = [all_E_by_band{b}, all_sca_by_band{b}, all_tau_by_band{b}, ...
                all_DOS_by_band{b}, all_V2_by_band{b}];
        dlmwrite(filename, data, 'delimiter', '\t', 'precision', '%.6e');
    end
end

fprintf('✓ All data saved to: %s\n', outputDir);
fprintf('\n=== Complete Processing Finished ===\n');