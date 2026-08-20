%% 1. Load the data files
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

fprintf('Files saved: %s_individual.dat and %s_average_filtered.dat\n', baseName, baseName);
