%% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';

%% Initialize for plotting - Create figures
figure1 = figure('Name', 'Scattering Rates');
fig1 = gcf;
fig1.Position(3:4) = [1200, 800];
hold on;

figure2 = figure('Name', 'Relaxation Times');
fig2 = gcf;
fig2.Position(3:4) = [1200, 800];
hold on;

numBands = 3;

% Initialize arrays for data collection - SEVEN METHODS
E_avg_all_M1 = []; sca_avg_all_M1 = []; tau_avg_all_M1 = [];
E_avg_all_M2 = []; sca_avg_all_M2 = []; tau_avg_all_M2 = [];
E_avg_all_M3 = []; sca_avg_all_M3 = []; tau_avg_all_M3 = [];
E_avg_all_M4 = []; sca_avg_all_M4 = []; tau_avg_all_M4 = [];
E_avg_all_M5 = []; sca_avg_all_M5 = []; tau_avg_all_M5 = [];
E_avg_all_M6 = []; sca_avg_all_M6 = []; tau_avg_all_M6 = [];
E_avg_all_M7 = []; sca_avg_all_M7 = []; tau_avg_all_M7 = [];

% Initialize per-band storage
E_avg_by_band_M1 = cell(numBands, 1); sca_avg_by_band_M1 = cell(numBands, 1); tau_avg_by_band_M1 = cell(numBands, 1);
E_avg_by_band_M2 = cell(numBands, 1); sca_avg_by_band_M2 = cell(numBands, 1); tau_avg_by_band_M2 = cell(numBands, 1);
E_avg_by_band_M3 = cell(numBands, 1); sca_avg_by_band_M3 = cell(numBands, 1); tau_avg_by_band_M3 = cell(numBands, 1);
E_avg_by_band_M4 = cell(numBands, 1); sca_avg_by_band_M4 = cell(numBands, 1); tau_avg_by_band_M4 = cell(numBands, 1);
E_avg_by_band_M5 = cell(numBands, 1); sca_avg_by_band_M5 = cell(numBands, 1); tau_avg_by_band_M5 = cell(numBands, 1);
E_avg_by_band_M6 = cell(numBands, 1); sca_avg_by_band_M6 = cell(numBands, 1); tau_avg_by_band_M6 = cell(numBands, 1);
E_avg_by_band_M7 = cell(numBands, 1); sca_avg_by_band_M7 = cell(numBands, 1); tau_avg_by_band_M7 = cell(numBands, 1);

% Store intermediate quantities
E_avg_by_band_V2g = cell(numBands, 1);
V2g_avg_by_band = cell(numBands, 1);
TDF_per_band_M6 = cell(numBands, 1);
TDF_IIS_per_band_M7 = cell(numBands, 1);

% Individual k-points
all_E_by_band = cell(numBands, 1);
all_sca_by_band = cell(numBands, 1);
all_tau_by_band = cell(numBands, 1);
all_DOS_by_band = cell(numBands, 1);
all_V2_by_band = cell(numBands, 1);

% Configuration
EF_index = 8;            % Index for EF=0 (8th column in dim 3)
T_index = 4;             % Temperature index (1st in dim 4 = 900K)

%% Check for required variables
fprintf('========================================\n');
fprintf('Checking required variables\n');
fprintf('========================================\n');

% Check TDF (total)
if exist('TDF', 'var') && isfield(TDF, 'xx')
    fprintf('✓ TDF.xx found [%d x %d x %d]\n', size(TDF.xx));
else
    error('TDF.xx not found!');
end

% Check TDF_ph (phonon total)
if exist('TDF_ph', 'var') && isfield(TDF_ph, 'xx')
    fprintf('✓ TDF_ph.xx found [%d x %d x %d]\n', size(TDF_ph.xx));
else
    error('TDF_ph.xx not found!');
end

% Check TDF_n (band-resolved total)
if exist('TDF_n', 'var') && isfield(TDF_n, 'xx')
    fprintf('✓ TDF_n.xx found [%d x %d x %d x %d]\n', size(TDF_n.xx));
else
    error('TDF_n.xx not found!');
end

% Check TDF_ph_n (band-resolved phonon)
if exist('TDF_ph_n', 'var') && isfield(TDF_ph_n, 'xx')
    fprintf('✓ TDF_ph_n.xx found [%d x %d x %d x %d]\n', size(TDF_ph_n.xx));
else
    error('TDF_ph_n.xx not found!');
end

fprintf('\n');

%% Extract TDF data at EF_index, T_index for M5 and M7
fprintf('========================================\n');
fprintf('Extracting TDF data at EF=%d, T=%d\n', EF_index, T_index);
fprintf('========================================\n');

% M5: Total IIS TDF from Matthiessen's rule
TDF_xx_total = TDF.xx(:, EF_index, T_index);
TDF_ph_xx_total = TDF_ph.xx(:, EF_index, T_index);

TDF_IIS_total = zeros(size(TDF_xx_total));
valid_mask_total = TDF_xx_total > 0 & TDF_ph_xx_total > 0 & TDF_xx_total < TDF_ph_xx_total;
TDF_IIS_total(valid_mask_total) = 1 ./ (1./TDF_xx_total(valid_mask_total) - 1./TDF_ph_xx_total(valid_mask_total));
fallback_mask_total = TDF_xx_total > 0 & ~valid_mask_total;
TDF_IIS_total(fallback_mask_total) = TDF_xx_total(fallback_mask_total);

fprintf('M5 - Total IIS TDF (Matthiessen):\n');
fprintf('  Total TDF:  [%.4e, %.4e]\n', min(TDF_xx_total(TDF_xx_total>0)), max(TDF_xx_total));
fprintf('  Phonon TDF: [%.4e, %.4e]\n', min(TDF_ph_xx_total(TDF_ph_xx_total>0)), max(TDF_ph_xx_total));
fprintf('  IIS TDF:    [%.4e, %.4e]\n', min(TDF_IIS_total(TDF_IIS_total>0)), max(TDF_IIS_total));
fprintf('  Matthiessen applied: %d/%d points\n\n', sum(valid_mask_total), length(valid_mask_total));

% M7: Band-resolved IIS TDF from Matthiessen's rule
TDF_n_bands = cell(numBands, 1);
TDF_ph_n_bands = cell(numBands, 1);
TDF_IIS_n_bands = cell(numBands, 1);

for b = 1:numBands
    TDF_n_bands{b} = squeeze(TDF_n.xx(:, b, EF_index, T_index));
    TDF_ph_n_bands{b} = squeeze(TDF_ph_n.xx(:, b, EF_index, T_index));
    
    TDF_IIS_b = zeros(size(TDF_n_bands{b}));
    valid_mask_b = TDF_n_bands{b} > 0 & TDF_ph_n_bands{b} > 0 & TDF_n_bands{b} < TDF_ph_n_bands{b};
    TDF_IIS_b(valid_mask_b) = 1 ./ (1./TDF_n_bands{b}(valid_mask_b) - 1./TDF_ph_n_bands{b}(valid_mask_b));
    fallback_mask_b = TDF_n_bands{b} > 0 & ~valid_mask_b;
    TDF_IIS_b(fallback_mask_b) = TDF_n_bands{b}(fallback_mask_b);
    
    TDF_IIS_n_bands{b} = TDF_IIS_b;
    
    pos_total = TDF_n_bands{b}(TDF_n_bands{b} > 0);
    pos_ph = TDF_ph_n_bands{b}(TDF_ph_n_bands{b} > 0);
    pos_IIS = TDF_IIS_b(TDF_IIS_b > 0);
    
    fprintf('M7 - Band %d IIS TDF (Matthiessen):\n', b);
    fprintf('  Total: [%.4e, %.4e], Phonon: [%.4e, %.4e]\n', min(pos_total), max(pos_total), min(pos_ph), max(pos_ph));
    fprintf('  IIS:   [%.4e, %.4e], Applied: %d/%d points\n\n', min(pos_IIS), max(pos_IIS), sum(valid_mask_b), length(valid_mask_b));
end

%% Get energy range from E_array
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

fprintf('\n========================================\n');
fprintf('Processing IIS data for EF=0, T=900K\n');
fprintf('Energy range: %.4f to %.4f eV\n', E_min, E_max);
fprintf('Number of energy points: %d\n', max_idx);
fprintf('========================================\n');
fprintf('\nSEVEN AVERAGING METHODS (all IIS-only for M4-M7):\n');
fprintf('M1: Simple <S>, tau = 1/<S>\n');
fprintf('M2: DOS-wt <S>, tau = 1/<S>\n');
fprintf('M3: DOS-wt <tau>, S = 1/<tau>\n');
fprintf('M4: Transport-wt <tau> (V²·g)\n');
fprintf('M5: tau = TDF_IIS_total / Σ(V²·g), TDF_IIS = [1/TDF - 1/TDF_ph]⁻¹\n');
fprintf('M6: tau = Σ(TDF_band_calc) / Σ(V²·g), TDF_band = Σ(τ_k·V²_k·g_k)\n');
fprintf('M7: tau = Σ(TDF_IIS_n) / Σ(V²·g), TDF_IIS_n = [1/TDF_n - 1/TDF_ph_n]⁻¹\n');

%% Main processing loop for M1-M4 and collecting for M5/M6/M7
for idx = 1:max_idx
    i = energy_indices(idx);
    
    E_row = [];
    sca_M1_row = []; tau_M1_row = [];
    sca_M2_row = []; tau_M2_row = [];
    sca_M3_row = []; tau_M3_row = [];
    sca_M4_row = []; tau_M4_row = [];
    DOS_row = [];
    
    TDF_band_at_E_M6 = zeros(1, numBands);
    
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
                    
                    % V²·g sum for this band
                    V2g_sum = sum(V2_k_positive .* DOS_k_positive);
                    
                    % M6: TDF_band = Σ(τ_k × V²_k × g_k)
                    TDF_band_val_M6 = sum(tau_k_positive .* V2_k_positive .* DOS_k_positive);
                    TDF_band_at_E_M6(b) = TDF_band_val_M6;
                    
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
                    
                    E_avg_by_band_V2g{b}(end+1) = E;
                    V2g_avg_by_band{b}(end+1) = V2g_sum;
                    
                    if isempty(TDF_per_band_M6{b})
                        TDF_per_band_M6{b} = [];
                    end
                    TDF_per_band_M6{b}(end+1) = TDF_band_val_M6;
                    
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
    
    % Overall M1-M4 and M6
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
        
        % M6: sum TDF across bands and divide by total V²g
        TDF_total_M6 = sum(TDF_band_at_E_M6);
        V2g_total = 0;
        for bb = 1:numBands
            if ~isempty(V2g_avg_by_band{bb})
                V2g_total = V2g_total + V2g_avg_by_band{bb}(end);
            end
        end
        
        if V2g_total > 0 && TDF_total_M6 > 0
            tau_M6_overall = TDF_total_M6 / V2g_total;
        else
            tau_M6_overall = NaN;
        end
        
        E_avg_all_M6(end+1) = E_avg;
        tau_avg_all_M6(end+1) = tau_M6_overall;
        sca_avg_all_M6(end+1) = 1 / tau_M6_overall;
    end
    
    if mod(idx, 50) == 0
        fprintf('  Processed %d/%d energy points...\n', idx, max_idx);
    end
end

N_overall = length(E_avg_all_M1);
fprintf('\nOverall energy points after processing: %d\n', N_overall);

%% METHOD 5: Total IIS TDF from Matthiessen's rule
fprintf('\n=== M5: Total IIS TDF = [1/TDF.xx - 1/TDF_ph.xx]^{-1} ===\n');

E_avg_all_M5 = zeros(N_overall, 1);
tau_avg_all_M5 = zeros(N_overall, 1);
sca_avg_all_M5 = zeros(N_overall, 1);

for j = 1:N_overall
    E_val = E_avg_all_M1(j);
    
    % Sum V²·g across bands
    V2g_total = 0;
    for b = 1:numBands
        if ~isempty(E_avg_by_band_V2g{b})
            [~, idx_b] = min(abs(E_avg_by_band_V2g{b} - E_val));
            if idx_b <= length(V2g_avg_by_band{b})
                V2g_total = V2g_total + V2g_avg_by_band{b}(idx_b);
            end
        end
    end
    
    % Interpolate TDF_IIS_total at this energy
    if exist('E_array', 'var') && length(E_array) == length(TDF_IIS_total)
        TDF_IIS_at_E = interp1(E_array, TDF_IIS_total, E_val, 'linear', 0);
    else
        [~, idx_e] = min(abs(E_array - E_val));
        TDF_IIS_at_E = TDF_IIS_total(min(idx_e, length(TDF_IIS_total)));
    end
    TDF_IIS_at_E = max(TDF_IIS_at_E, 0);
    
    if V2g_total > 0 && TDF_IIS_at_E > 0
        tau_M5 = TDF_IIS_at_E / V2g_total;
    else
        tau_M5 = NaN;
    end
    
    E_avg_all_M5(j) = E_val;
    tau_avg_all_M5(j) = tau_M5;
    sca_avg_all_M5(j) = 1 / tau_M5;
end

valid_M5 = ~isnan(tau_avg_all_M5) & tau_avg_all_M5 > 0;
E_avg_all_M5_valid = E_avg_all_M5(valid_M5);
tau_avg_all_M5_valid = tau_avg_all_M5(valid_M5);
sca_avg_all_M5_valid = sca_avg_all_M5(valid_M5);

fprintf('M5 complete. %d valid points out of %d.\n', sum(valid_M5), N_overall);

%% METHOD 7: Band-resolved IIS TDF from Matthiessen's rule
fprintf('\n=== M7: Band-resolved IIS TDF = [1/TDF_n.xx - 1/TDF_ph_n.xx]^{-1} ===\n');
fprintf('    tau = Σ(TDF_IIS_n) / Σ(V²·g)\n');

E_avg_all_M7 = zeros(N_overall, 1);
tau_avg_all_M7 = zeros(N_overall, 1);
sca_avg_all_M7 = zeros(N_overall, 1);

% Initialize per-band storage for M7
for b = 1:numBands
    TDF_IIS_per_band_M7{b} = zeros(N_overall, 1);
end

for j = 1:N_overall
    E_val = E_avg_all_M1(j);
    
    % Sum V²·g across bands
    V2g_total = 0;
    for b = 1:numBands
        if ~isempty(E_avg_by_band_V2g{b})
            [~, idx_b] = min(abs(E_avg_by_band_V2g{b} - E_val));
            if idx_b <= length(V2g_avg_by_band{b})
                V2g_total = V2g_total + V2g_avg_by_band{b}(idx_b);
            end
        end
    end
    
    % Sum band-resolved IIS TDFs
    TDF_IIS_sum = 0;
    for b = 1:numBands
        if exist('E_array', 'var') && length(E_array) == length(TDF_IIS_n_bands{b})
            TDF_IIS_val = interp1(E_array, TDF_IIS_n_bands{b}, E_val, 'linear', 0);
        else
            [~, idx_e] = min(abs(E_array - E_val));
            if idx_e <= length(TDF_IIS_n_bands{b})
                TDF_IIS_val = TDF_IIS_n_bands{b}(idx_e);
            else
                TDF_IIS_val = 0;
            end
        end
        TDF_IIS_val = max(TDF_IIS_val, 0);
        
        TDF_IIS_per_band_M7{b}(j) = TDF_IIS_val;
        TDF_IIS_sum = TDF_IIS_sum + TDF_IIS_val;
    end
    
    if V2g_total > 0 && TDF_IIS_sum > 0
        tau_M7 = TDF_IIS_sum / V2g_total;
    else
        tau_M7 = NaN;
    end
    
    E_avg_all_M7(j) = E_val;
    tau_avg_all_M7(j) = tau_M7;
    sca_avg_all_M7(j) = 1 / tau_M7;
end

valid_M7 = ~isnan(tau_avg_all_M7) & tau_avg_all_M7 > 0;
E_avg_all_M7_valid = E_avg_all_M7(valid_M7);
tau_avg_all_M7_valid = tau_avg_all_M7(valid_M7);
sca_avg_all_M7_valid = sca_avg_all_M7(valid_M7);

fprintf('M7 complete. %d valid points out of %d.\n', sum(valid_M7), N_overall);

% Clean M6
valid_M6 = ~isnan(tau_avg_all_M6) & tau_avg_all_M6 > 0;
E_avg_all_M6_valid = E_avg_all_M6(valid_M6);
tau_avg_all_M6_valid = tau_avg_all_M6(valid_M6);
sca_avg_all_M6_valid = sca_avg_all_M6(valid_M6);

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
        plot(E_avg_all_M5_valid, sca_avg_all_M5_valid, 'm-', 'LineWidth', 2.5, 'DisplayName', 'M5: 1/<tau>_{IIS-total}');
    end
    if exist('E_avg_all_M6_valid', 'var')
        plot(E_avg_all_M6_valid, sca_avg_all_M6_valid, 'c-', 'LineWidth', 2.5, 'DisplayName', 'M6: 1/<tau>_{IIS-calc}');
    end
    if exist('E_avg_all_M7_valid', 'var')
        plot(E_avg_all_M7_valid, sca_avg_all_M7_valid, 'Color', [0.8 0.4 0], 'LineWidth', 3, ...
             'DisplayName', 'M7: 1/<tau>_{IIS-bands}');
    end
end

set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([E_min, E_max]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Scattering Rate [fs^{-1}]', 'FontSize', 12);
title('ZrNiSn - Scattering Rates (7 Methods)', 'FontSize', 14);
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
        plot(E_avg_all_M5_valid, tau_avg_all_M5_valid, 'm-', 'LineWidth', 2.5, 'DisplayName', 'M5: TDF_{IIS-total}/V²g');
    end
    if exist('E_avg_all_M6_valid', 'var')
        plot(E_avg_all_M6_valid, tau_avg_all_M6_valid, 'c-', 'LineWidth', 2.5, 'DisplayName', 'M6: ΣTDF_{calc}/V²g');
    end
    if exist('E_avg_all_M7_valid', 'var')
        plot(E_avg_all_M7_valid, tau_avg_all_M7_valid, 'Color', [0.8 0.4 0], 'LineWidth', 3, ...
             'DisplayName', 'M7: ΣTDF_{IIS-bands}/V²g');
    end
end

set(gca, 'YScale', 'log');
ylim([1e-3, 1e1]);
xlim([E_min, E_max]);
grid on;
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Relaxation Time \tau [fs]', 'FontSize', 12);
title('ZrNiSn - Relaxation Times (7 Methods)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 8);

%% STATISTICS
fprintf('\n========================================\n');
fprintf('COMPREHENSIVE METHOD COMPARISON\n');
fprintf('========================================\n');

fprintf('\n--- Mean tau Values ---\n');
fprintf('M1: %.4f fs\n', mean(tau_avg_all_M1(~isnan(tau_avg_all_M1) & tau_avg_all_M1 > 0)));
fprintf('M2: %.4f fs\n', mean(tau_avg_all_M2(~isnan(tau_avg_all_M2) & tau_avg_all_M2 > 0)));
fprintf('M3: %.4f fs\n', mean(tau_avg_all_M3(~isnan(tau_avg_all_M3) & tau_avg_all_M3 > 0)));
fprintf('M4: %.4f fs\n', mean(tau_avg_all_M4(~isnan(tau_avg_all_M4) & tau_avg_all_M4 > 0)));
fprintf('M5: %.4f fs\n', mean(tau_avg_all_M5_valid));
fprintf('M6: %.4f fs\n', mean(tau_avg_all_M6_valid));
fprintf('M7: %.4f fs\n', mean(tau_avg_all_M7_valid));

% Interpolate all to M7 grid for comparison
fprintf('\n--- Ratios Relative to M7 ---\n');
tau_M1_interp = interp1(E_avg_all_M1, tau_avg_all_M1, E_avg_all_M7_valid, 'linear');
tau_M2_interp = interp1(E_avg_all_M2, tau_avg_all_M2, E_avg_all_M7_valid, 'linear');
tau_M3_interp = interp1(E_avg_all_M3, tau_avg_all_M3, E_avg_all_M7_valid, 'linear');
tau_M4_interp = interp1(E_avg_all_M4, tau_avg_all_M4, E_avg_all_M7_valid, 'linear');
tau_M5_interp = interp1(E_avg_all_M5_valid, tau_avg_all_M5_valid, E_avg_all_M7_valid, 'linear');
tau_M6_interp = interp1(E_avg_all_M6_valid, tau_avg_all_M6_valid, E_avg_all_M7_valid, 'linear');

fprintf('M1/M7: %.4f\n', nanmean(tau_M1_interp ./ tau_avg_all_M7_valid));
fprintf('M2/M7: %.4f\n', nanmean(tau_M2_interp ./ tau_avg_all_M7_valid));
fprintf('M3/M7: %.4f\n', nanmean(tau_M3_interp ./ tau_avg_all_M7_valid));
fprintf('M4/M7: %.4f\n', nanmean(tau_M4_interp ./ tau_avg_all_M7_valid));
fprintf('M5/M7: %.4f\n', nanmean(tau_M5_interp ./ tau_avg_all_M7_valid));
fprintf('M6/M7: %.4f\n', nanmean(tau_M6_interp ./ tau_avg_all_M7_valid));

% Also compare M4, M5, M6 among themselves
fprintf('\n--- IIS Method Cross-Comparison ---\n');
fprintf('M4/M5: %.4f\n', nanmean(tau_M4_interp ./ tau_M5_interp));
fprintf('M4/M6: %.4f\n', nanmean(tau_M4_interp ./ tau_M6_interp));
fprintf('M5/M6: %.4f\n', nanmean(tau_M5_interp ./ tau_M6_interp));

%% SAVE DATA
fprintf('\n=== Saving data ===\n');
baseName = 'ZrNiSn_IIS_900K_7methods';
outputDir = 'IIS_xmgrace_7methods_900K';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save overall 7-method comparison
filename = fullfile(outputDir, sprintf('%s_overall.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# 7-method comparison: tau values at EF=0, T=900K\n');
fprintf(fid, '# All IIS-only for M4-M7\n');
fprintf(fid, '# M4: Transport-wt <tau> = Sigma(tau*V^2*g) / Sigma(V^2*g)\n');
fprintf(fid, '# M5: TDF_IIS_total / Sigma(V^2*g), TDF_IIS = [1/TDF - 1/TDF_ph]^{-1}\n');
fprintf(fid, '# M6: Sigma(TDF_band_calc) / Sigma(V^2*g), TDF_band = Sigma(tau_k*V^2_k*g_k)\n');
fprintf(fid, '# M7: Sigma(TDF_IIS_n) / Sigma(V^2*g), TDF_IIS_n = [1/TDF_n - 1/TDF_ph_n]^{-1}\n');
fprintf(fid, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5  Tau_M6  Tau_M7\n');

n_pts = length(E_avg_all_M7_valid);
for j = 1:n_pts
    E_val = E_avg_all_M7_valid(j);
    tau1 = interp1(E_avg_all_M1, tau_avg_all_M1, E_val, 'linear', NaN);
    tau2 = interp1(E_avg_all_M2, tau_avg_all_M2, E_val, 'linear', NaN);
    tau3 = interp1(E_avg_all_M3, tau_avg_all_M3, E_val, 'linear', NaN);
    tau4 = interp1(E_avg_all_M4, tau_avg_all_M4, E_val, 'linear', NaN);
    tau5 = interp1(E_avg_all_M5_valid, tau_avg_all_M5_valid, E_val, 'linear', NaN);
    tau6 = interp1(E_avg_all_M6_valid, tau_avg_all_M6_valid, E_val, 'linear', NaN);
    tau7 = tau_avg_all_M7_valid(j);
    
    fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
            E_val, tau1, tau2, tau3, tau4, tau5, tau6, tau7);
end
fclose(fid);

% Save band-resolved IIS TDF data for M7
filename_M7 = fullfile(outputDir, sprintf('%s_band_IIS_TDF_M7.dat', baseName));
fid_M7 = fopen(filename_M7, 'w');
fprintf(fid_M7, '# Band-resolved IIS TDF for M7\n');
fprintf(fid_M7, '# TDF_IIS_n = [1/TDF_n.xx - 1/TDF_ph_n.xx]^{-1} at EF=0, T=900K\n');
fprintf(fid_M7, '# Energy  TDF_IIS_B1  TDF_IIS_B2  TDF_IIS_B3  TDF_IIS_Sum  V2g_Total  Tau_M7\n');

for j = 1:n_pts
    E_val = E_avg_all_M7_valid(j);
    
    % Get V2g total
    V2g_total = 0;
    for b = 1:numBands
        if ~isempty(E_avg_by_band_V2g{b})
            [~, idx_b] = min(abs(E_avg_by_band_V2g{b} - E_val));
            if idx_b <= length(V2g_avg_by_band{b})
                V2g_total = V2g_total + V2g_avg_by_band{b}(idx_b);
            end
        end
    end
    
    fprintf(fid_M7, '%.6f', E_val);
    TDF_sum = 0;
    for b = 1:numBands
        tdf_val = TDF_IIS_per_band_M7{b}(j);
        fprintf(fid_M7, '  %.6e', tdf_val);
        TDF_sum = TDF_sum + tdf_val;
    end
    fprintf(fid_M7, '  %.6e  %.6e  %.6e\n', TDF_sum, V2g_total, tau_avg_all_M7_valid(j));
end
fclose(fid_M7);

% Save M5 IIS TDF data
filename_M5 = fullfile(outputDir, sprintf('%s_total_IIS_TDF_M5.dat', baseName));
fid_M5 = fopen(filename_M5, 'w');
fprintf(fid_M5, '# Total IIS TDF for M5\n');
fprintf(fid_M5, '# TDF_IIS = [1/TDF.xx - 1/TDF_ph.xx]^{-1} at EF=0, T=900K\n');
fprintf(fid_M5, '# Energy  TDF_total  TDF_ph  TDF_IIS  V2g_Total  Tau_M5\n');

for j = 1:length(E_avg_all_M5_valid)
    E_val = E_avg_all_M5_valid(j);
    
    % Interpolate TDF values
    if exist('E_array', 'var')
        TDF_total_at_E = interp1(E_array, TDF_xx_total, E_val, 'linear', 0);
        TDF_ph_at_E = interp1(E_array, TDF_ph_xx_total, E_val, 'linear', 0);
        TDF_IIS_at_E = interp1(E_array, TDF_IIS_total, E_val, 'linear', 0);
    else
        TDF_total_at_E = 0;
        TDF_ph_at_E = 0;
        TDF_IIS_at_E = 0;
    end
    
    % Get V2g total
    V2g_total = 0;
    for b = 1:numBands
        if ~isempty(E_avg_by_band_V2g{b})
            [~, idx_b] = min(abs(E_avg_by_band_V2g{b} - E_val));
            if idx_b <= length(V2g_avg_by_band{b})
                V2g_total = V2g_total + V2g_avg_by_band{b}(idx_b);
            end
        end
    end
    
    fprintf(fid_M5, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
            E_val, TDF_total_at_E, TDF_ph_at_E, TDF_IIS_at_E, V2g_total, tau_avg_all_M5_valid(j));
end
fclose(fid_M5);

fprintf('✓ Data saved to: %s\n', outputDir);
fprintf('  - %s_overall.dat (all 7 methods)\n', baseName);
fprintf('  - %s_band_IIS_TDF_M7.dat (per-band IIS TDF)\n', baseName);
fprintf('  - %s_total_IIS_TDF_M5.dat (total IIS TDF with components)\n', baseName);
fprintf('\n=== Complete ===\n');
