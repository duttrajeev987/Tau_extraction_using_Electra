% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_matthiessen_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_screenedPOP_ZrNiSn_kScan_holes.mat';
%load 'RelaxTimes_POP_ZrNiSn_kScan_holes.mat'
%load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';
% Initialize for plotting
figure;
fig = gcf;
fig.Position(3:4) = [550, 400];
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
        if ~isempty(taus_POP(i,b).x) && size(taus_POP(i,b).x, 3) >= 1
            % Get energy for this point
            E = state_ID(i,b).E;
            tau_3D = taus_POP(i,b).x;
            
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
    dlmwrite('scattering_individual_ZrNiSn_POP_300K.dat', individual_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('Saved %d individual k-points to scattering_individual_ZrNiSn_POP_300K.dat\n', length(all_E));
else
    fprintf('Warning: No individual k-point data to save\n');
end

% 2. Save energy-averaged data (for line plots
if ~isempty(E_avg_all)
    avg_data = [E_avg_all(:), sca_avg_all(:)];
    dlmwrite('scattering_average_ZrNiSn_POP_300K.dat', avg_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('Saved %d averaged points to scattering_average_ZrNiSn_POP_300K.dat\n', length(E_avg_all));
else
    fprintf('Warning: No averaged data to save\n');
end

% 3. Save with band information (3 columns: E, sca, band)
if ~isempty(all_E)
    band_data = [all_E, all_sca, all_band];
    dlmwrite('scattering_with_bands.dat', band_data, ...
        'delimiter', '\t', 'precision', '%.6e');
    fprintf('Saved band-labeled data to scattering_with_bands_ZrNiSn_POP_300K.dat\n');
end

% 4. Save with headers (xmgrace ignores lines starting with #)
if ~isempty(E_avg_all)
    fid = fopen('scattering_formatted_ZrNiSn_POP_300K.dat', 'w');
    fprintf(fid, '# Energy (eV)\tScattering Rate (fs^-1)\n');
    fprintf(fid, '# Data extracted at EF=0, T=300K\n');
    fprintf(fid, '# Filter: Only positive tau values included\n');
    for i = 1:length(E_avg_all)
        fprintf(fid, '%.6f\t%.6e\n', E_avg_all(i), sca_avg_all(i));
    end
    fclose(fid);
    fprintf('Saved formatted data to scattering_formatted_ZrNiSn_POP_300K.dat\n');
end

% 5. Save MATLAB workspace for further analysis
if ~isempty(E_avg_all) || ~isempty(all_E)
    save('scattering_analysis_ZrNiSn_POP_300.mat', 'E_avg_all', 'sca_avg_all', ...
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
