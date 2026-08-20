% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_matthiessen_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_screenedPOP_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_POP_ZrNiSn_kScan_holes.mat'
load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';
% Initialize for plotting
figure;
fig = gcf;
fig.Position(3:4) = [550, 400];
hold on;

numBands = 3;
palette = [
    0.1216 0.46%% 1. Load the data files
load 'TE_ZrNiSn_kScan_holes.mat'; 
load 'RelaxTimes_matthiessen_ZrNiSn_kScan_holes.mat'; 

%% 2. Configuration & Parameters
targetTempIdx = 1;      % 1=300K, 2=400K, 3=500K, 4=600K, etc.
threshold = 1e15;       % Threshold: Ignore rates > 10^15 for the average
tempVal = 300 + (targetTempIdx-1)*100;

% Visualization Setup
figure;
fig = gcf;
fig.Position(3:4) = [600, 450];
hold on;
palette = [0.1216 0.4667 0.7059]; % Professional Blue

% Initialize data containers
E_avg_all   = [];          
sca_avg_all = [];
all_E       = [];              
all_sca     = [];            

%% 3. Main Processing Loop
% Iterating through all energy slices (standard index 1 to 451)
for i = 1:length(taus_matth)
    
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
            
            % Filter 2: Apply the 10^15 threshold for the AVERAGE calculation
            % We keep the "noisy" points for the background scatter but exclude from mean
            avg_mask = sca_k <= threshold;
            sca_k_for_avg = sca_k(avg_mask);
            
            % Collect ALL individual points for plotting/saving
            all_E = [all_E; repmat(E, length(sca_k), 1)];
            all_sca = [all_sca; sca_k];
            
            % Calculate and store average only if valid points exist below threshold
            if ~isempty(sca_k_for_avg)
                sca_slice_avg = mean(sca_k_for_avg);
                
                E_avg_all(end+1)   = E;
                sca_avg_all(end+1) = sca_slice_avg;
                
                % Plot individual k-points (light/faded)
                scatter(repmat(E, size(sca_k)), sca_k, 20, palette, 'filled', 'MarkerFaceAlpha', 0.1);
                
                % Plot the filtered average point (bold diamond)
                scatter(E, sca_slice_avg, 45, 'r', 'filled', 'd'); 
            end
        end
    end
end

%% 4. Finalize Plot
avg_line = plot(E_avg_all, sca_avg_all, 'r-', 'LineWidth', 2.5);

set(gca, 'YScale', 'log');
grid on;
ylim([1e11, 2e15]); % Adjusted slightly to see the threshold limit
xlim([0, 0.45]);
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Scattering Rate \tau^{-1} [s^{-1}]', 'FontSize', 12);
title(['ZrNiSn Matthiessen | T = ', num2str(tempVal), 'K | Filter: < 10^{15}'], 'FontSize', 14);

legend([avg_line], {['Average Rate (Filtered < ', num2str(threshold, '%.0e'), ')']}, 'Location', 'best');

%% 5. Export Data for xmgrace
fprintf('\n=== Exporting Matthiessen Data (%dK) ===\n', tempVal);
baseName = sprintf('ZrNiSn_Matth_%dK', tempVal);

% A. Save Individual K-points
if ~isempty(all_E)
    dlmwrite([baseName, '_individual.dat'], [all_E, all_sca], 'delimiter', '\t', 'precision', '%.6e');
end

% B. Save Filtered Averages (The Trend Line)
if ~isempty(E_avg_all)
    dlmwrite([baseName, '_average_filtered.dat'], [E_avg_all(:), sca_avg_all(:)], 'delimiter', '\t', 'precision', '%.6e');
end

% C. Formatted file for readability
fid = fopen([baseName, '_summary.dat'], 'w');
fprintf(fid, '# ZrNiSn Scattering Data (Matthiessen Rule)\n');
fprintf(fid, '# Temperature: %d K\n', tempVal);
fprintf(fid, '# Filter: Scattering rates > %.1e excluded from average\n', threshold);
fprintf(fid, '# Energy (eV)\tAvg_Scat_Rate (s^-1)\n');
for i = 1:length(E_avg_all)
    fprintf(fid, '%.6f\t%.6e\n', E_avg_all(i), sca_avg_all(i));
end
fclose(fid);

fprintf('Files saved: %s_individual.dat and %s_average_filtered.dat\n', baseName, baseName);67 0.7059;  % blue
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
all_band = [];           % For individual k-points

% Main processing loop
for i = 1:451
    E_row = [];
    sca_row = [];
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x, 3) >= 1
            % Get energy for this point
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D, 1) >= 8 && size(tau_3D, 2) >= 1
                % Extract tau(EF=0, T=300K, all k-points)
                tau_k = squeeze(tau_3D(8, 1, :));
                %tau_k = tau_3D(:, 1);% for IISiessen values 1 stands for temperature at 1st place in T_array
                
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
                    if i == 1 && ~h_created(b)
                        h(b) = htmp;
                        h_created(b) = true;
                    end
                    
                    % Collect for total average over bands
                    E_row(end+1) = E;
                    sca_row(end+1) = sca_band_avg;
                else
                    fprintf("No positive tau values for i=%d, b=%d\n", i, b);
                end
            end
        end
    end
    
    % After looping over bands, take the overall average
    if ~isempty(E_row)
        E_avg_all(end+1) = mean(E_row);
        sca_avg_all(end+1) = mean(sca_row);
    end
end

% Overlay average line
avg_line = plot(E_avg_all, sca_avg_all, 'k-', 'LineWidth', 2);

% Format the plot
set(gca, 'YScale', 'log');
ylim([1e11, 1e15]);
xlim([0, 0.45]);

xlabel('Energy [eV]');
ylabel('\tau^{-1} [fs^{-1}]');
title('ZrNiSn - EF=0, T=300K, all k-points');

% Create legend - FIXED VERSION
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
    lgd.FontSize = 14;
    lgd.Location = 'best';
else
    warning('No valid data to create legend');
end

%% Save data for xmgrace plotting
fprintf('\n=== Saving data for xmgrace ===\n');

% 1. Save all individual k-points (for scatter plots)
if ~isempty(all_E)
    individual_data = [all_E, all_sca];
    dlmwrite('scattering_individual_ZrNiSn_IIS_300K.dat', individual_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('Saved %d individual k-points to scattering_individual_ZrNiSn_IIS_300K.dat\n', length(all_E));
else
    fprintf('Warning: No individual k-point data to save\n');
end

% 2. Save energy-averaged data (for line plots
if ~isempty(E_avg_all)
    avg_data = [E_avg_all(:), sca_avg_all(:)];
    dlmwrite('scattering_average_ZrNiSn_IIS_300K.dat', avg_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('Saved %d averaged points to scattering_average_ZrNiSn_IIS_300K.dat\n', length(E_avg_all));
else
    fprintf('Warning: No averaged data to save\n');
end

% 3. Save with band information (3 columns: E, sca, band)
if ~isempty(all_E)
    band_data = [all_E, all_sca, all_band];
    dlmwrite('scattering_with_bands.dat', band_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('Saved band-labeled data to scattering_with_bands_ZrNiSn_IIS_300K.dat\n');
end

% 4. Save with headers (xmgrace ignores lines starting with #)
if ~isempty(E_avg_all)
    fid = fopen('scattering_formatted_ZrNiSn_IIS_300K.dat', 'w');
    fprintf(fid, '# Energy (eV)\tScattering Rate (fs^-1)\n');
    fprintf(fid, '# Data extracted at EF=0, T=300K\n');
    fprintf(fid, '# Filter: Only positive tau values included\n');
    for i = 1:length(E_avg_all)
        fprintf(fid, '%.6f\t%.6e\n', E_avg_all(i), sca_avg_all(i));
    end
    fclose(fid);
    fprintf('Saved formatted data to scattering_formatted_ZrNiSn_IIS_300K.dat\n');
end

% 5. Save MATLAB workspace for further analysis
if ~isempty(E_avg_all) || ~isempty(all_E)
    save('scattering_analysis_ZrNiSn_IIS_300.mat', 'E_avg_all', 'sca_avg_all', ...
        'all_E', 'all_sca', 'all_band', 'h_created');
    fprintf('Saved workspace to scattering_analysis.mat\n');
end

fprintf('\n=== Summary ===\n');
fprintf('Processed %d energy points\n', 451);
fprintf('Found valid data in %d bands\n', sum(h_created));
if ~isempty(all_E)
    fprintf('Total individual k-points: %d\n', length(all_E));
end
if ~isempty(E_avg_all)
    fprintf('Average curve points: %d\n', length(E_avg_all));
end