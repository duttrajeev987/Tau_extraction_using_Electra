%% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_separate_ZrNiSn_kScan_holes.mat';

%% Configuration
numBands = 3;

% Define temperature configuration
T_indices = [1, 2, 3, 4, 5, 6, 7];  % Temperature indices (1=300K, 7=900K)
T_labels = {'300K', '400K', '500K', '600K', '700K', '800K', '900K'};

% Mechanism indices
ADP_idx = 1;  % First dimension index for ADP
ODP_idx = 2;  % First dimension index for ODP

% Check maximum temperature index available
if exist('taus', 'var')
    sample_idx = find(~cellfun(@isempty, {taus(:,1).x}), 1);
    if ~isempty(sample_idx)
        max_T_available = size(taus(sample_idx,1).x, 3);
        fprintf('Maximum temperature index available: %d\n', max_T_available);
        T_indices = T_indices(T_indices <= max_T_available);
        T_labels = T_labels(1:length(T_indices));
    end
else
    error('taus not found! Check if the variable exists in the loaded file.');
end

fprintf('Processing %d temperatures: %s\n', length(T_indices), strjoin(T_labels, ', '));
fprintf('ADP mechanism index: %d, ODP mechanism index: %d\n', ADP_idx, ODP_idx);

%% Check for required variables
fprintf('========================================\n');
fprintf('Checking required variables\n');
fprintf('========================================\n');

% Check taus structure
if exist('taus', 'var')
    fprintf('✓ taus found [%d x %d]\n', size(taus, 1), size(taus, 2));
    for i = 1:min(3, size(taus, 1))
        for b = 1:min(numBands, size(taus, 2))
            if ~isempty(taus(i,b).x)
                fprintf('  taus(%d,%d).x size: [', i, b);
                fprintf('%d ', size(taus(i,b).x));
                fprintf(']\n');
                for mech = 1:2
                    sample_tau = squeeze(taus(i,b).x(mech, :, 1));
                    sample_tau = sample_tau(sample_tau > 0);
                    if ~isempty(sample_tau)
                        mech_name = {'ADP', 'ODP'};
                        fprintf('    %s sample tau: [%.2e, %.2e]\n', mech_name{mech}, min(sample_tau), max(sample_tau));
                    end
                end
                break;
            end
        end
    end
else
    error('taus not found!');
end

% Check TDF_sep structure
TDF_available = false;
has_band_resolved_ADP = false;
has_band_resolved_ODP = false;

if exist('TDF_sep', 'var')
    fprintf('✓ TDF_sep found\n');
    
    if isfield(TDF_sep, 'ADP') && isfield(TDF_sep.ADP, 'xx')
        fprintf('✓ TDF_sep.ADP.xx found [%d x %d x %d]\n', size(TDF_sep.ADP.xx));
        has_ADP = true;
    else
        warning('TDF_sep.ADP.xx not found!');
        has_ADP = false;
    end
    
    if isfield(TDF_sep, 'ODP') && isfield(TDF_sep.ODP, 'xx')
        fprintf('✓ TDF_sep.ODP.xx found [%d x %d x %d]\n', size(TDF_sep.ODP.xx));
        has_ODP = true;
    else
        warning('TDF_sep.ODP.xx not found!');
        has_ODP = false;
    end
    
    TDF_available = has_ADP && has_ODP;
    
    if isfield(TDF_sep, 'ADP_n') && isfield(TDF_sep.ADP_n, 'xx')
        has_band_resolved_ADP = true;
        fprintf('✓ TDF_sep.ADP_n.xx found\n');
    end
    if isfield(TDF_sep, 'ODP_n') && isfield(TDF_sep.ODP_n, 'xx')
        has_band_resolved_ODP = true;
        fprintf('✓ TDF_sep.ODP_n.xx found\n');
    end
    
    has_band_resolved = has_band_resolved_ADP && has_band_resolved_ODP;
end

fprintf('\n');

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
    warning('E_array not found.');
end

%% Main loop over temperatures
for temp_idx = 1:length(T_indices)
    T_index = T_indices(temp_idx);
    T_label = T_labels{temp_idx};
    
    fprintf('\n========================================\n');
    fprintf('PROCESSING TEMPERATURE: %s (Index: %d)\n', T_label, T_index);
    fprintf('========================================\n');
    
    %% Initialize figures
    figure1 = figure('Name', sprintf('Scattering Rates - ADP+ODP+Matthiessen - %s', T_label));
    fig1 = gcf;
    fig1.Position(3:4) = [1200, 800];
    hold on;
    
    figure2 = figure('Name', sprintf('Relaxation Times - ADP+ODP+Matthiessen - %s', T_label));
    fig2 = gcf;
    fig2.Position(3:4) = [1200, 800];
    hold on;
    
    % ============ ADP Arrays ============
    E_avg_all_ADP = []; 
    sca_avg_all_M1_ADP = []; tau_avg_all_M1_ADP = [];
    sca_avg_all_M2_ADP = []; tau_avg_all_M2_ADP = [];
    sca_avg_all_M3_ADP = []; tau_avg_all_M3_ADP = [];
    sca_avg_all_M4_ADP = []; tau_avg_all_M4_ADP = [];
    sca_avg_all_M5_ADP = []; tau_avg_all_M5_ADP = [];
    sca_avg_all_M6_ADP = []; tau_avg_all_M6_ADP = [];
    
    E_avg_by_band_ADP = cell(numBands, 1);
    V2g_avg_by_band_ADP = cell(numBands, 1);
    TDF_per_band_M6_ADP = cell(numBands, 1);
    
    % ============ ODP Arrays ============
    E_avg_all_ODP = []; 
    sca_avg_all_M1_ODP = []; tau_avg_all_M1_ODP = [];
    sca_avg_all_M2_ODP = []; tau_avg_all_M2_ODP = [];
    sca_avg_all_M3_ODP = []; tau_avg_all_M3_ODP = [];
    sca_avg_all_M4_ODP = []; tau_avg_all_M4_ODP = [];
    sca_avg_all_M5_ODP = []; tau_avg_all_M5_ODP = [];
    sca_avg_all_M6_ODP = []; tau_avg_all_M6_ODP = [];
    
    E_avg_by_band_ODP = cell(numBands, 1);
    V2g_avg_by_band_ODP = cell(numBands, 1);
    TDF_per_band_M6_ODP = cell(numBands, 1);
    
    % ============ Matthiessen Arrays ============
    E_avg_all_MATTH = []; 
    sca_avg_all_M1_MATTH = []; tau_avg_all_M1_MATTH = [];
    sca_avg_all_M2_MATTH = []; tau_avg_all_M2_MATTH = [];
    sca_avg_all_M3_MATTH = []; tau_avg_all_M3_MATTH = [];
    sca_avg_all_M4_MATTH = []; tau_avg_all_M4_MATTH = [];
    sca_avg_all_M5_MATTH = []; tau_avg_all_M5_MATTH = [];
    sca_avg_all_M6_MATTH = []; tau_avg_all_M6_MATTH = [];
    sca_avg_all_M7_MATTH = []; tau_avg_all_M7_MATTH = [];
    
    E_avg_by_band_MATTH = cell(numBands, 1);
    V2g_avg_by_band_MATTH = cell(numBands, 1);
    TDF_per_band_M6_MATTH = cell(numBands, 1);
    TDF_matth_per_band_M7 = cell(numBands, 1);
    
    % Initialize per-band cells
    for b = 1:numBands
        E_avg_by_band_ADP{b} = [];
        V2g_avg_by_band_ADP{b} = [];
        TDF_per_band_M6_ADP{b} = [];
        E_avg_by_band_ODP{b} = [];
        V2g_avg_by_band_ODP{b} = [];
        TDF_per_band_M6_ODP{b} = [];
        E_avg_by_band_MATTH{b} = [];
        V2g_avg_by_band_MATTH{b} = [];
        TDF_per_band_M6_MATTH{b} = [];
    end
    
    %% Extract TDF data
    TDF_ADP_total = [];
    TDF_ODP_total = [];
    TDF_matth_total = [];
    TDF_ADP_n_bands = cell(numBands, 1);
    TDF_ODP_n_bands = cell(numBands, 1);
    TDF_matth_n_bands = cell(numBands, 1);
    
    if TDF_available
        ADP_data = TDF_sep.ADP.xx;
        ODP_data = TDF_sep.ODP.xx;
        
        % Extract ADP TDF
        if ndims(ADP_data) == 3
            for test_ef = 1:size(ADP_data, 2)
                if size(ADP_data, 3) >= T_index
                    test_data = ADP_data(:, test_ef, T_index);
                    if any(test_data > 0)
                        TDF_ADP_total = test_data;
                        break;
                    end
                end
            end
        elseif ndims(ADP_data) == 2
            if size(ADP_data, 2) >= T_index
                TDF_ADP_total = ADP_data(:, T_index);
            end
        end
        
        % Extract ODP TDF
        if ndims(ODP_data) == 3
            for test_ef = 1:size(ODP_data, 2)
                if size(ODP_data, 3) >= T_index
                    test_data = ODP_data(:, test_ef, T_index);
                    if any(test_data > 0)
                        TDF_ODP_total = test_data;
                        break;
                    end
                end
            end
        elseif ndims(ODP_data) == 2
            if size(ODP_data, 2) >= T_index
                TDF_ODP_total = ODP_data(:, T_index);
            end
        end
        
        % Calculate Matthiessen TDF
        if ~isempty(TDF_ADP_total) && ~isempty(TDF_ODP_total) && ...
           any(TDF_ADP_total > 0) && any(TDF_ODP_total > 0)
            
            valid_mask = TDF_ADP_total > 0 & TDF_ODP_total > 0;
            TDF_matth_total = zeros(size(TDF_ADP_total));
            TDF_matth_total(valid_mask) = 1 ./ (1./TDF_ADP_total(valid_mask) + 1./TDF_ODP_total(valid_mask));
            
            fprintf('TDF data extracted for %s:\n', T_label);
            fprintf('  ADP: [%.4e, %.4e], Non-zero: %d\n', ...
                    min(TDF_ADP_total(TDF_ADP_total>0)), max(TDF_ADP_total), sum(TDF_ADP_total>0));
            fprintf('  ODP: [%.4e, %.4e], Non-zero: %d\n', ...
                    min(TDF_ODP_total(TDF_ODP_total>0)), max(TDF_ODP_total), sum(TDF_ODP_total>0));
            fprintf('  Matth: [%.4e, %.4e], Non-zero: %d\n', ...
                    min(TDF_matth_total(TDF_matth_total>0)), max(TDF_matth_total), sum(TDF_matth_total>0));
        end
        
        % Band-resolved TDF (simplified - skip if causing issues)
        if has_band_resolved
            try
                if has_band_resolved_ADP
                    ADP_n_data = TDF_sep.ADP_n.xx;
                    if ndims(ADP_n_data) == 4
                        for b = 1:numBands
                            if size(ADP_n_data, 4) >= T_index
                                for test_ef = 1:size(ADP_n_data, 3)
                                    test_data = squeeze(ADP_n_data(:, b, test_ef, T_index));
                                    if any(test_data > 0)
                                        TDF_ADP_n_bands{b} = test_data;
                                        break;
                                    end
                                end
                            end
                        end
                    end
                end
                
                if has_band_resolved_ODP
                    ODP_n_data = TDF_sep.ODP_n.xx;
                    if ndims(ODP_n_data) == 4
                        for b = 1:numBands
                            if size(ODP_n_data, 4) >= T_index
                                for test_ef = 1:size(ODP_n_data, 3)
                                    test_data = squeeze(ODP_n_data(:, b, test_ef, T_index));
                                    if any(test_data > 0)
                                        TDF_ODP_n_bands{b} = test_data;
                                        break;
                                    end
                                end
                            end
                        end
                    end
                end
                
                for b = 1:numBands
                    if ~isempty(TDF_ADP_n_bands{b}) && ~isempty(TDF_ODP_n_bands{b})
                        valid_mask = TDF_ADP_n_bands{b} > 0 & TDF_ODP_n_bands{b} > 0;
                        TDF_matth_n_bands{b} = zeros(size(TDF_ADP_n_bands{b}));
                        TDF_matth_n_bands{b}(valid_mask) = 1 ./ (1./TDF_ADP_n_bands{b}(valid_mask) + 1./TDF_ODP_n_bands{b}(valid_mask));
                    end
                end
            catch
                warning('Band-resolved TDF extraction failed, M7 will use total TDF.');
            end
        end
    end
    
    %% Main processing loop
    fprintf('Processing energy points...\n');
    
    tau_is_seconds = false;
    
    for idx = 1:max_idx
        i = energy_indices(idx);
        
        % Rows for each mechanism
        E_row_ADP = []; DOS_row_ADP = [];
        sca_M1_ADP_row = []; tau_M1_ADP_row = [];
        sca_M2_ADP_row = []; tau_M2_ADP_row = [];
        sca_M3_ADP_row = []; tau_M3_ADP_row = [];
        sca_M4_ADP_row = []; tau_M4_ADP_row = [];
        TDF_band_at_E_M6_ADP = zeros(1, numBands);
        
        E_row_ODP = []; DOS_row_ODP = [];
        sca_M1_ODP_row = []; tau_M1_ODP_row = [];
        sca_M2_ODP_row = []; tau_M2_ODP_row = [];
        sca_M3_ODP_row = []; tau_M3_ODP_row = [];
        sca_M4_ODP_row = []; tau_M4_ODP_row = [];
        TDF_band_at_E_M6_ODP = zeros(1, numBands);
        
        E_row_MATTH = []; DOS_row_MATTH = [];
        sca_M1_MATTH_row = []; tau_M1_MATTH_row = [];
        sca_M2_MATTH_row = []; tau_M2_MATTH_row = [];
        sca_M3_MATTH_row = []; tau_M3_MATTH_row = [];
        sca_M4_MATTH_row = []; tau_M4_MATTH_row = [];
        TDF_band_at_E_M6_MATTH = zeros(1, numBands);
        
        for b = 1:numBands
            if size(taus, 1) >= i && size(taus, 2) >= b
                if ~isempty(taus(i,b).x) && size(taus(i,b).x, 3) >= T_index
                    E = state_ID(i,b).E;
                    tau_3D = taus(i,b).x;
                    
                    % Extract ADP and ODP tau
                    tau_ADP_k = squeeze(tau_3D(ADP_idx, :, T_index));
                    tau_ODP_k = squeeze(tau_3D(ODP_idx, :, T_index));
                    
                    % Ensure column vectors
                    tau_ADP_k = tau_ADP_k(:);
                    tau_ODP_k = tau_ODP_k(:);
                    
                    % Get number of k-points
                    nk = length(tau_ADP_k);
                    
                    % Calculate Matthiessen combined tau
                    valid_both = tau_ADP_k > 0 & tau_ODP_k > 0;
                    tau_MATTH_k = zeros(nk, 1);
                    tau_MATTH_k(valid_both) = 1 ./ (1./tau_ADP_k(valid_both) + 1./tau_ODP_k(valid_both));
                    
                    % Check if tau is in seconds
                    if ~tau_is_seconds && idx == 1 && b == 1
                        max_tau = max([tau_ADP_k(tau_ADP_k > 0); tau_ODP_k(tau_ODP_k > 0)]);
                        if ~isempty(max_tau) && max_tau < 1e-12 && max_tau > 0
                            tau_is_seconds = true;
                            fprintf('  NOTE: Tau values in seconds, converting to fs (×1e15)\n');
                        end
                    end
                    
                    if tau_is_seconds
                        tau_ADP_k = tau_ADP_k * 1e15;
                        tau_ODP_k = tau_ODP_k * 1e15;
                        tau_MATTH_k = tau_MATTH_k * 1e15;
                    end
                    
                    % Extract velocity - ensure correct size
                    if isfield(state_ID(i,b), 'V') && ~isempty(state_ID(i,b).V)
                        V_k = state_ID(i,b).V(:);
                        if length(V_k) == nk
                            V2_k = V_k.^2;
                        else
                            V2_k = ones(nk, 1);
                        end
                    else
                        V2_k = ones(nk, 1);
                    end
                    
                    % Extract DOS - ensure correct size
                    if isfield(state_ID(i,b), 'DOS') && ~isempty(state_ID(i,b).DOS)
                        DOS_k = state_ID(i,b).DOS(:);
                        if length(DOS_k) ~= nk
                            DOS_k = ones(nk, 1);
                        end
                    else
                        DOS_k = ones(nk, 1);
                    end
                    
                    % ============ Process ADP ============
                    positive_mask_ADP = tau_ADP_k > 0 & V2_k > 0 & DOS_k > 0;
                    if any(positive_mask_ADP)
                        tau_ADP_pos = tau_ADP_k(positive_mask_ADP);
                        V2_ADP_pos = V2_k(positive_mask_ADP);
                        DOS_ADP_pos = DOS_k(positive_mask_ADP);
                        sca_ADP_k = 1 ./ tau_ADP_pos;
                        
                        % M1-M4
                        sca_M1_ADP = mean(sca_ADP_k);
                        tau_M1_ADP = 1 / sca_M1_ADP;
                        
                        sca_M2_ADP = sum(sca_ADP_k .* DOS_ADP_pos) / sum(DOS_ADP_pos);
                        tau_M2_ADP = 1 / sca_M2_ADP;
                        
                        tau_M3_ADP = sum(tau_ADP_pos .* DOS_ADP_pos) / sum(DOS_ADP_pos);
                        sca_M3_ADP = 1 / tau_M3_ADP;
                        
                        transport_wt_ADP = V2_ADP_pos .* DOS_ADP_pos;
                        tau_M4_ADP = sum(tau_ADP_pos .* transport_wt_ADP) / sum(transport_wt_ADP);
                        sca_M4_ADP = 1 / tau_M4_ADP;
                        
                        V2g_sum_ADP = sum(transport_wt_ADP);
                        TDF_band_val_M6_ADP = sum(tau_ADP_pos .* transport_wt_ADP);
                        TDF_band_at_E_M6_ADP(b) = TDF_band_val_M6_ADP;
                        
                        % Store per-band
                        E_avg_by_band_ADP{b}(end+1) = E;
                        V2g_avg_by_band_ADP{b}(end+1) = V2g_sum_ADP;
                        TDF_per_band_M6_ADP{b}(end+1) = TDF_band_val_M6_ADP;
                        
                        % Collect rows
                        E_row_ADP(end+1) = E;
                        sca_M1_ADP_row(end+1) = sca_M1_ADP; tau_M1_ADP_row(end+1) = tau_M1_ADP;
                        sca_M2_ADP_row(end+1) = sca_M2_ADP; tau_M2_ADP_row(end+1) = tau_M2_ADP;
                        sca_M3_ADP_row(end+1) = sca_M3_ADP; tau_M3_ADP_row(end+1) = tau_M3_ADP;
                        sca_M4_ADP_row(end+1) = sca_M4_ADP; tau_M4_ADP_row(end+1) = tau_M4_ADP;
                        DOS_row_ADP(end+1) = sum(DOS_ADP_pos);
                    end
                    
                    % ============ Process ODP ============
                    positive_mask_ODP = tau_ODP_k > 0 & V2_k > 0 & DOS_k > 0;
                    if any(positive_mask_ODP)
                        tau_ODP_pos = tau_ODP_k(positive_mask_ODP);
                        V2_ODP_pos = V2_k(positive_mask_ODP);
                        DOS_ODP_pos = DOS_k(positive_mask_ODP);
                        sca_ODP_k = 1 ./ tau_ODP_pos;
                        
                        sca_M1_ODP = mean(sca_ODP_k);
                        tau_M1_ODP = 1 / sca_M1_ODP;
                        
                        sca_M2_ODP = sum(sca_ODP_k .* DOS_ODP_pos) / sum(DOS_ODP_pos);
                        tau_M2_ODP = 1 / sca_M2_ODP;
                        
                        tau_M3_ODP = sum(tau_ODP_pos .* DOS_ODP_pos) / sum(DOS_ODP_pos);
                        sca_M3_ODP = 1 / tau_M3_ODP;
                        
                        transport_wt_ODP = V2_ODP_pos .* DOS_ODP_pos;
                        tau_M4_ODP = sum(tau_ODP_pos .* transport_wt_ODP) / sum(transport_wt_ODP);
                        sca_M4_ODP = 1 / tau_M4_ODP;
                        
                        V2g_sum_ODP = sum(transport_wt_ODP);
                        TDF_band_val_M6_ODP = sum(tau_ODP_pos .* transport_wt_ODP);
                        TDF_band_at_E_M6_ODP(b) = TDF_band_val_M6_ODP;
                        
                        E_avg_by_band_ODP{b}(end+1) = E;
                        V2g_avg_by_band_ODP{b}(end+1) = V2g_sum_ODP;
                        TDF_per_band_M6_ODP{b}(end+1) = TDF_band_val_M6_ODP;
                        
                        E_row_ODP(end+1) = E;
                        sca_M1_ODP_row(end+1) = sca_M1_ODP; tau_M1_ODP_row(end+1) = tau_M1_ODP;
                        sca_M2_ODP_row(end+1) = sca_M2_ODP; tau_M2_ODP_row(end+1) = tau_M2_ODP;
                        sca_M3_ODP_row(end+1) = sca_M3_ODP; tau_M3_ODP_row(end+1) = tau_M3_ODP;
                        sca_M4_ODP_row(end+1) = sca_M4_ODP; tau_M4_ODP_row(end+1) = tau_M4_ODP;
                        DOS_row_ODP(end+1) = sum(DOS_ODP_pos);
                    end
                    
                    % ============ Process Matthiessen ============
                    positive_mask_MATTH = tau_MATTH_k > 0 & V2_k > 0 & DOS_k > 0;
                    if any(positive_mask_MATTH)
                        tau_MATTH_pos = tau_MATTH_k(positive_mask_MATTH);
                        V2_MATTH_pos = V2_k(positive_mask_MATTH);
                        DOS_MATTH_pos = DOS_k(positive_mask_MATTH);
                        sca_MATTH_k = 1 ./ tau_MATTH_pos;
                        
                        sca_M1_MATTH = mean(sca_MATTH_k);
                        tau_M1_MATTH = 1 / sca_M1_MATTH;
                        
                        sca_M2_MATTH = sum(sca_MATTH_k .* DOS_MATTH_pos) / sum(DOS_MATTH_pos);
                        tau_M2_MATTH = 1 / sca_M2_MATTH;
                        
                        tau_M3_MATTH = sum(tau_MATTH_pos .* DOS_MATTH_pos) / sum(DOS_MATTH_pos);
                        sca_M3_MATTH = 1 / tau_M3_MATTH;
                        
                        transport_wt_MATTH = V2_MATTH_pos .* DOS_MATTH_pos;
                        tau_M4_MATTH = sum(tau_MATTH_pos .* transport_wt_MATTH) / sum(transport_wt_MATTH);
                        sca_M4_MATTH = 1 / tau_M4_MATTH;
                        
                        V2g_sum_MATTH = sum(transport_wt_MATTH);
                        TDF_band_val_M6_MATTH = sum(tau_MATTH_pos .* transport_wt_MATTH);
                        TDF_band_at_E_M6_MATTH(b) = TDF_band_val_M6_MATTH;
                        
                        E_avg_by_band_MATTH{b}(end+1) = E;
                        V2g_avg_by_band_MATTH{b}(end+1) = V2g_sum_MATTH;
                        TDF_per_band_M6_MATTH{b}(end+1) = TDF_band_val_M6_MATTH;
                        
                        E_row_MATTH(end+1) = E;
                        sca_M1_MATTH_row(end+1) = sca_M1_MATTH; tau_M1_MATTH_row(end+1) = tau_M1_MATTH;
                        sca_M2_MATTH_row(end+1) = sca_M2_MATTH; tau_M2_MATTH_row(end+1) = tau_M2_MATTH;
                        sca_M3_MATTH_row(end+1) = sca_M3_MATTH; tau_M3_MATTH_row(end+1) = tau_M3_MATTH;
                        sca_M4_MATTH_row(end+1) = sca_M4_MATTH; tau_M4_MATTH_row(end+1) = tau_M4_MATTH;
                        DOS_row_MATTH(end+1) = sum(DOS_MATTH_pos);
                    end
                end
            end
        end
        
        % ============ Overall ADP averages ============
        if ~isempty(E_row_ADP) && sum(DOS_row_ADP) > 0
            E_avg = sum(E_row_ADP .* DOS_row_ADP) / sum(DOS_row_ADP);
            
            E_avg_all_ADP(end+1) = E_avg;
            sca_avg_all_M1_ADP(end+1) = sum(sca_M1_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            tau_avg_all_M1_ADP(end+1) = sum(tau_M1_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            sca_avg_all_M2_ADP(end+1) = sum(sca_M2_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            tau_avg_all_M2_ADP(end+1) = sum(tau_M2_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            sca_avg_all_M3_ADP(end+1) = sum(sca_M3_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            tau_avg_all_M3_ADP(end+1) = sum(tau_M3_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            sca_avg_all_M4_ADP(end+1) = sum(sca_M4_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            tau_avg_all_M4_ADP(end+1) = sum(tau_M4_ADP_row .* DOS_row_ADP) / sum(DOS_row_ADP);
            
            % M6 ADP
            TDF_total_M6_ADP = sum(TDF_band_at_E_M6_ADP);
            V2g_total_ADP = 0;
            for bb = 1:numBands
                if ~isempty(V2g_avg_by_band_ADP{bb})
                    V2g_total_ADP = V2g_total_ADP + V2g_avg_by_band_ADP{bb}(end);
                end
            end
            if V2g_total_ADP > 0 && TDF_total_M6_ADP > 0
                tau_M6_ADP = TDF_total_M6_ADP / V2g_total_ADP;
            else
                tau_M6_ADP = NaN;
            end
            sca_avg_all_M6_ADP(end+1) = 1 / tau_M6_ADP;
            tau_avg_all_M6_ADP(end+1) = tau_M6_ADP;
        end
        
        % ============ Overall ODP averages ============
        if ~isempty(E_row_ODP) && sum(DOS_row_ODP) > 0
            E_avg = sum(E_row_ODP .* DOS_row_ODP) / sum(DOS_row_ODP);
            
            E_avg_all_ODP(end+1) = E_avg;
            sca_avg_all_M1_ODP(end+1) = sum(sca_M1_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            tau_avg_all_M1_ODP(end+1) = sum(tau_M1_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            sca_avg_all_M2_ODP(end+1) = sum(sca_M2_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            tau_avg_all_M2_ODP(end+1) = sum(tau_M2_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            sca_avg_all_M3_ODP(end+1) = sum(sca_M3_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            tau_avg_all_M3_ODP(end+1) = sum(tau_M3_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            sca_avg_all_M4_ODP(end+1) = sum(sca_M4_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            tau_avg_all_M4_ODP(end+1) = sum(tau_M4_ODP_row .* DOS_row_ODP) / sum(DOS_row_ODP);
            
            TDF_total_M6_ODP = sum(TDF_band_at_E_M6_ODP);
            V2g_total_ODP = 0;
            for bb = 1:numBands
                if ~isempty(V2g_avg_by_band_ODP{bb})
                    V2g_total_ODP = V2g_total_ODP + V2g_avg_by_band_ODP{bb}(end);
                end
            end
            if V2g_total_ODP > 0 && TDF_total_M6_ODP > 0
                tau_M6_ODP = TDF_total_M6_ODP / V2g_total_ODP;
            else
                tau_M6_ODP = NaN;
            end
            sca_avg_all_M6_ODP(end+1) = 1 / tau_M6_ODP;
            tau_avg_all_M6_ODP(end+1) = tau_M6_ODP;
        end
        
        % ============ Overall Matthiessen averages ============
        if ~isempty(E_row_MATTH) && sum(DOS_row_MATTH) > 0
            E_avg = sum(E_row_MATTH .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            
            E_avg_all_MATTH(end+1) = E_avg;
            sca_avg_all_M1_MATTH(end+1) = sum(sca_M1_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            tau_avg_all_M1_MATTH(end+1) = sum(tau_M1_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            sca_avg_all_M2_MATTH(end+1) = sum(sca_M2_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            tau_avg_all_M2_MATTH(end+1) = sum(tau_M2_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            sca_avg_all_M3_MATTH(end+1) = sum(sca_M3_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            tau_avg_all_M3_MATTH(end+1) = sum(tau_M3_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            sca_avg_all_M4_MATTH(end+1) = sum(sca_M4_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            tau_avg_all_M4_MATTH(end+1) = sum(tau_M4_MATTH_row .* DOS_row_MATTH) / sum(DOS_row_MATTH);
            
            TDF_total_M6_MATTH = sum(TDF_band_at_E_M6_MATTH);
            V2g_total_MATTH = 0;
            for bb = 1:numBands
                if ~isempty(V2g_avg_by_band_MATTH{bb})
                    V2g_total_MATTH = V2g_total_MATTH + V2g_avg_by_band_MATTH{bb}(end);
                end
            end
            if V2g_total_MATTH > 0 && TDF_total_M6_MATTH > 0
                tau_M6_MATTH = TDF_total_M6_MATTH / V2g_total_MATTH;
            else
                tau_M6_MATTH = NaN;
            end
            sca_avg_all_M6_MATTH(end+1) = 1 / tau_M6_MATTH;
            tau_avg_all_M6_MATTH(end+1) = tau_M6_MATTH;
        end
        
        if mod(idx, 50) == 0
            fprintf('  Processed %d/%d energy points...\n', idx, max_idx);
        end
    end
    
    N_ADP = length(E_avg_all_ADP);
    N_ODP = length(E_avg_all_ODP);
    N_MATTH = length(E_avg_all_MATTH);
    fprintf('\nEnergy points: ADP=%d, ODP=%d, Matthiessen=%d\n', N_ADP, N_ODP, N_MATTH);
    
    %% M5: TDF-based calculations
    
    % ADP M5
    if ~isempty(TDF_ADP_total) && any(TDF_ADP_total > 0) && N_ADP > 0
        sca_avg_all_M5_ADP = zeros(N_ADP, 1);
        tau_avg_all_M5_ADP = zeros(N_ADP, 1);
        for j = 1:N_ADP
            E_val = E_avg_all_ADP(j);
            V2g_total = 0;
            for b = 1:numBands
                if ~isempty(E_avg_by_band_ADP{b})
                    [~, idx_b] = min(abs(E_avg_by_band_ADP{b} - E_val));
                    if idx_b <= length(V2g_avg_by_band_ADP{b})
                        V2g_total = V2g_total + V2g_avg_by_band_ADP{b}(idx_b);
                    end
                end
            end
            if exist('E_array', 'var') && length(E_array) == length(TDF_ADP_total)
                TDF_at_E = interp1(E_array, TDF_ADP_total, E_val, 'linear', 0);
            else
                [~, idx_e] = min(abs(E_array - E_val));
                TDF_at_E = TDF_ADP_total(min(idx_e, length(TDF_ADP_total)));
            end
            TDF_at_E = max(TDF_at_E, 0);
            if tau_is_seconds, TDF_at_E = TDF_at_E * 1e15; end
            
            if V2g_total > 0 && TDF_at_E > 0
                tau_avg_all_M5_ADP(j) = TDF_at_E / V2g_total;
                sca_avg_all_M5_ADP(j) = 1 / tau_avg_all_M5_ADP(j);
            else
                tau_avg_all_M5_ADP(j) = NaN;
                sca_avg_all_M5_ADP(j) = NaN;
            end
        end
    end
    
    % ODP M5
    if ~isempty(TDF_ODP_total) && any(TDF_ODP_total > 0) && N_ODP > 0
        sca_avg_all_M5_ODP = zeros(N_ODP, 1);
        tau_avg_all_M5_ODP = zeros(N_ODP, 1);
        for j = 1:N_ODP
            E_val = E_avg_all_ODP(j);
            V2g_total = 0;
            for b = 1:numBands
                if ~isempty(E_avg_by_band_ODP{b})
                    [~, idx_b] = min(abs(E_avg_by_band_ODP{b} - E_val));
                    if idx_b <= length(V2g_avg_by_band_ODP{b})
                        V2g_total = V2g_total + V2g_avg_by_band_ODP{b}(idx_b);
                    end
                end
            end
            if exist('E_array', 'var') && length(E_array) == length(TDF_ODP_total)
                TDF_at_E = interp1(E_array, TDF_ODP_total, E_val, 'linear', 0);
            else
                [~, idx_e] = min(abs(E_array - E_val));
                TDF_at_E = TDF_ODP_total(min(idx_e, length(TDF_ODP_total)));
            end
            TDF_at_E = max(TDF_at_E, 0);
            if tau_is_seconds, TDF_at_E = TDF_at_E * 1e15; end
            
            if V2g_total > 0 && TDF_at_E > 0
                tau_avg_all_M5_ODP(j) = TDF_at_E / V2g_total;
                sca_avg_all_M5_ODP(j) = 1 / tau_avg_all_M5_ODP(j);
            else
                tau_avg_all_M5_ODP(j) = NaN;
                sca_avg_all_M5_ODP(j) = NaN;
            end
        end
    end
    
    % Matthiessen M5 & M7
    if ~isempty(TDF_matth_total) && any(TDF_matth_total > 0) && N_MATTH > 0
        sca_avg_all_M5_MATTH = zeros(N_MATTH, 1);
        tau_avg_all_M5_MATTH = zeros(N_MATTH, 1);
        sca_avg_all_M7_MATTH = zeros(N_MATTH, 1);
        tau_avg_all_M7_MATTH = zeros(N_MATTH, 1);
        
        for b = 1:numBands
            TDF_matth_per_band_M7{b} = zeros(N_MATTH, 1);
        end
        
        for j = 1:N_MATTH
            E_val = E_avg_all_MATTH(j);
            
            V2g_total = 0;
            for b = 1:numBands
                if ~isempty(E_avg_by_band_MATTH{b})
                    [~, idx_b] = min(abs(E_avg_by_band_MATTH{b} - E_val));
                    if idx_b <= length(V2g_avg_by_band_MATTH{b})
                        V2g_total = V2g_total + V2g_avg_by_band_MATTH{b}(idx_b);
                    end
                end
            end
            
            % M5
            if exist('E_array', 'var') && length(E_array) == length(TDF_matth_total)
                TDF_at_E = interp1(E_array, TDF_matth_total, E_val, 'linear', 0);
            else
                [~, idx_e] = min(abs(E_array - E_val));
                TDF_at_E = TDF_matth_total(min(idx_e, length(TDF_matth_total)));
            end
            TDF_at_E = max(TDF_at_E, 0);
            if tau_is_seconds, TDF_at_E = TDF_at_E * 1e15; end
            
            if V2g_total > 0 && TDF_at_E > 0
                tau_avg_all_M5_MATTH(j) = TDF_at_E / V2g_total;
                sca_avg_all_M5_MATTH(j) = 1 / tau_avg_all_M5_MATTH(j);
            else
                tau_avg_all_M5_MATTH(j) = NaN;
                sca_avg_all_M5_MATTH(j) = NaN;
            end
            
            % M7
            TDF_sum = 0;
            has_M7_data = false;
            for b = 1:numBands
                if ~isempty(TDF_matth_n_bands{b}) && any(TDF_matth_n_bands{b} > 0)
                    has_M7_data = true;
                    if exist('E_array', 'var') && length(E_array) == length(TDF_matth_n_bands{b})
                        TDF_val = interp1(E_array, TDF_matth_n_bands{b}, E_val, 'linear', 0);
                    else
                        [~, idx_e] = min(abs(E_array - E_val));
                        if idx_e <= length(TDF_matth_n_bands{b})
                            TDF_val = TDF_matth_n_bands{b}(idx_e);
                        else
                            TDF_val = 0;
                        end
                    end
                    TDF_val = max(TDF_val, 0);
                elseif ~isempty(TDF_matth_total) && any(TDF_matth_total > 0)
                    if exist('E_array', 'var')
                        TDF_val = interp1(E_array, TDF_matth_total, E_val, 'linear', 0) / numBands;
                    else
                        [~, idx_e] = min(abs(E_array - E_val));
                        TDF_val = TDF_matth_total(min(idx_e, length(TDF_matth_total))) / numBands;
                    end
                else
                    TDF_val = 0;
                end
                if tau_is_seconds, TDF_val = TDF_val * 1e15; end
                
                TDF_matth_per_band_M7{b}(j) = TDF_val;
                TDF_sum = TDF_sum + TDF_val;
            end
            
            if V2g_total > 0 && TDF_sum > 0
                tau_avg_all_M7_MATTH(j) = TDF_sum / V2g_total;
                sca_avg_all_M7_MATTH(j) = 1 / tau_avg_all_M7_MATTH(j);
            else
                tau_avg_all_M7_MATTH(j) = NaN;
                sca_avg_all_M7_MATTH(j) = NaN;
            end
        end
    end
    
    %% FIGURES (simplified)
    figure(figure1);
    if ~isempty(E_avg_all_ADP)
        plot(E_avg_all_ADP, sca_avg_all_M1_ADP, 'b-', 'LineWidth', 1.5, 'DisplayName', 'ADP M1');
        plot(E_avg_all_ADP, sca_avg_all_M4_ADP, 'b--', 'LineWidth', 1.5, 'DisplayName', 'ADP M4');
        plot(E_avg_all_ODP, sca_avg_all_M1_ODP, 'r-', 'LineWidth', 1.5, 'DisplayName', 'ODP M1');
        plot(E_avg_all_ODP, sca_avg_all_M4_ODP, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ODP M4');
        plot(E_avg_all_MATTH, sca_avg_all_M4_MATTH, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Matth M4');
        if exist('sca_avg_all_M5_MATTH', 'var')
            plot(E_avg_all_MATTH, sca_avg_all_M5_MATTH, 'm-', 'LineWidth', 2, 'DisplayName', 'Matth M5');
        end
        if exist('sca_avg_all_M7_MATTH', 'var')
            plot(E_avg_all_MATTH, sca_avg_all_M7_MATTH, 'Color', [0.8 0.4 0], 'LineWidth', 3, 'DisplayName', 'Matth M7');
        end
    end
    set(gca, 'YScale', 'log'); ylim auto; xlim([E_min, E_max]); grid on;
    xlabel('Energy [eV]'); ylabel('Scattering Rate [fs^{-1}]');
    title(sprintf('ZrNiSn - Scattering Rates (%s)', T_label));
    legend('Location', 'best', 'FontSize', 7);
    
    figure(figure2);
    if ~isempty(E_avg_all_ADP)
        plot(E_avg_all_ADP, tau_avg_all_M1_ADP, 'b-', 'LineWidth', 1.5, 'DisplayName', 'ADP M1');
        plot(E_avg_all_ADP, tau_avg_all_M4_ADP, 'b--', 'LineWidth', 1.5, 'DisplayName', 'ADP M4');
        plot(E_avg_all_ODP, tau_avg_all_M1_ODP, 'r-', 'LineWidth', 1.5, 'DisplayName', 'ODP M1');
        plot(E_avg_all_ODP, tau_avg_all_M4_ODP, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ODP M4');
        plot(E_avg_all_MATTH, tau_avg_all_M4_MATTH, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Matth M4');
        if exist('tau_avg_all_M5_MATTH', 'var')
            plot(E_avg_all_MATTH, tau_avg_all_M5_MATTH, 'm-', 'LineWidth', 2, 'DisplayName', 'Matth M5');
        end
        if exist('tau_avg_all_M7_MATTH', 'var')
            plot(E_avg_all_MATTH, tau_avg_all_M7_MATTH, 'Color', [0.8 0.4 0], 'LineWidth', 3, 'DisplayName', 'Matth M7');
        end
    end
    set(gca, 'YScale', 'log'); ylim auto; xlim([E_min, E_max]); grid on;
    xlabel('Energy [eV]'); ylabel('Relaxation Time \tau [fs]');
    title(sprintf('ZrNiSn - Relaxation Times (%s)', T_label));
    legend('Location', 'best', 'FontSize', 7);
    
    %% STATISTICS
    fprintf('\n=== STATISTICS FOR %s ===\n', T_label);
    fprintf('\n--- ADP ---\n');
    fprintf('M1: %.4e, M2: %.4e, M3: %.4e, M4: %.4e, M6: %.4e fs\n', ...
        nanmean(tau_avg_all_M1_ADP), nanmean(tau_avg_all_M2_ADP), nanmean(tau_avg_all_M3_ADP), ...
        nanmean(tau_avg_all_M4_ADP), nanmean(tau_avg_all_M6_ADP));
    fprintf('\n--- ODP ---\n');
    fprintf('M1: %.4e, M2: %.4e, M3: %.4e, M4: %.4e, M6: %.4e fs\n', ...
        nanmean(tau_avg_all_M1_ODP), nanmean(tau_avg_all_M2_ODP), nanmean(tau_avg_all_M3_ODP), ...
        nanmean(tau_avg_all_M4_ODP), nanmean(tau_avg_all_M6_ODP));
    fprintf('\n--- Matthiessen ---\n');
    fprintf('M1: %.4e, M2: %.4e, M3: %.4e, M4: %.4e, M6: %.4e', ...
        nanmean(tau_avg_all_M1_MATTH), nanmean(tau_avg_all_M2_MATTH), nanmean(tau_avg_all_M3_MATTH), ...
        nanmean(tau_avg_all_M4_MATTH), nanmean(tau_avg_all_M6_MATTH));
    if exist('tau_avg_all_M5_MATTH', 'var')
        fprintf(', M5: %.4e', nanmean(tau_avg_all_M5_MATTH));
    end
    if exist('tau_avg_all_M7_MATTH', 'var')
        fprintf(', M7: %.4e', nanmean(tau_avg_all_M7_MATTH));
    end
    fprintf('\n');
    
    %% SAVE DATA - Three separate directories
    fprintf('\n=== Saving data for %s ===\n', T_label);
    
    baseDir = sprintf('separate_mechanisms_%s', T_label);
    dirADP = fullfile(baseDir, 'ADP');
    dirODP = fullfile(baseDir, 'ODP');
    dirMATTH = fullfile(baseDir, 'Matthiessen');
    
    if ~exist(dirADP, 'dir'), mkdir(dirADP); end
    if ~exist(dirODP, 'dir'), mkdir(dirODP); end
    if ~exist(dirMATTH, 'dir'), mkdir(dirMATTH); end
    
    % Save ADP
    filename_ADP = fullfile(dirADP, sprintf('ZrNiSn_ADP_%s.dat', T_label));
    fid_ADP = fopen(filename_ADP, 'w');
    fprintf(fid_ADP, '# ADP relaxation times for ZrNiSn at %s\n', T_label);
    fprintf(fid_ADP, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5  Tau_M6\n');
    for j = 1:N_ADP
        tau5 = NaN;
        if exist('tau_avg_all_M5_ADP', 'var') && j <= length(tau_avg_all_M5_ADP)
            tau5 = tau_avg_all_M5_ADP(j);
        end
        fprintf(fid_ADP, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                E_avg_all_ADP(j), tau_avg_all_M1_ADP(j), tau_avg_all_M2_ADP(j), ...
                tau_avg_all_M3_ADP(j), tau_avg_all_M4_ADP(j), tau5, tau_avg_all_M6_ADP(j));
    end
    fclose(fid_ADP);
    fprintf('  ADP -> %s\n', dirADP);
    
    % Save ODP
    filename_ODP = fullfile(dirODP, sprintf('ZrNiSn_ODP_%s.dat', T_label));
    fid_ODP = fopen(filename_ODP, 'w');
    fprintf(fid_ODP, '# ODP relaxation times for ZrNiSn at %s\n', T_label);
    fprintf(fid_ODP, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5  Tau_M6\n');
    for j = 1:N_ODP
        tau5 = NaN;
        if exist('tau_avg_all_M5_ODP', 'var') && j <= length(tau_avg_all_M5_ODP)
            tau5 = tau_avg_all_M5_ODP(j);
        end
        fprintf(fid_ODP, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                E_avg_all_ODP(j), tau_avg_all_M1_ODP(j), tau_avg_all_M2_ODP(j), ...
                tau_avg_all_M3_ODP(j), tau_avg_all_M4_ODP(j), tau5, tau_avg_all_M6_ODP(j));
    end
    fclose(fid_ODP);
    fprintf('  ODP -> %s\n', dirODP);
    
    % Save Matthiessen
    filename_MATTH = fullfile(dirMATTH, sprintf('ZrNiSn_Matthiessen_%s.dat', T_label));
    fid_MATTH = fopen(filename_MATTH, 'w');
    fprintf(fid_MATTH, '# Matthiessen (ADP+ODP) relaxation times for ZrNiSn at %s\n', T_label);
    fprintf(fid_MATTH, '# 1/tau_matth = 1/tau_ADP + 1/tau_ODP\n');
    fprintf(fid_MATTH, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5  Tau_M6  Tau_M7\n');
    for j = 1:N_MATTH
        tau5 = NaN; tau7 = NaN;
        if exist('tau_avg_all_M5_MATTH', 'var') && j <= length(tau_avg_all_M5_MATTH)
            tau5 = tau_avg_all_M5_MATTH(j);
        end
        if exist('tau_avg_all_M7_MATTH', 'var') && j <= length(tau_avg_all_M7_MATTH)
            tau7 = tau_avg_all_M7_MATTH(j);
        end
        fprintf(fid_MATTH, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                E_avg_all_MATTH(j), tau_avg_all_M1_MATTH(j), tau_avg_all_M2_MATTH(j), ...
                tau_avg_all_M3_MATTH(j), tau_avg_all_M4_MATTH(j), tau5, ...
                tau_avg_all_M6_MATTH(j), tau7);
    end
    fclose(fid_MATTH);
    fprintf('  Matthiessen -> %s\n', dirMATTH);
    
end

fprintf('\n========================================\n');
fprintf('ALL TEMPERATURES PROCESSED SUCCESSFULLY\n');
fprintf('========================================\n');
fprintf('Data saved in separate directories for ADP, ODP, and Matthiessen\n');
fprintf('\n=== Complete ===\n');