%% SAVE DATA
fprintf('\n=== Saving data ===\n');
baseName = 'ZrNiSn_IIS_300K_TDF_corrected';
outputDir = 'IIS_xmgrace_TDF_corrected';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% SAVE OVERALL DATA
filename = fullfile(outputDir, sprintf('%s_overall_5methods.dat', baseName));
fid = fopen(filename, 'w');
fprintf(fid, '# Complete 5-method comparison - Overall (all bands combined)\n');
fprintf(fid, '# M5: tau = TDF.xx(E, EF=0, T=300K) / (V^2(E)*g(E))\n');
fprintf(fid, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5  Sca_M1  Sca_M2  Sca_M3  Sca_M4  Sca_M5\n');

% Use M5 energy grid if available
if exist('E_avg_all_M5_valid', 'var')
    n_pts = length(E_avg_all_M5_valid);
    for j = 1:n_pts
        E_val = E_avg_all_M5_valid(j);
        tau1 = interp1(E_avg_all_M1, tau_avg_all_M1, E_val, 'linear', NaN);
        tau2 = interp1(E_avg_all_M2, tau_avg_all_M2, E_val, 'linear', NaN);
        tau3 = interp1(E_avg_all_M3, tau_avg_all_M3, E_val, 'linear', NaN);
        tau4 = interp1(E_avg_all_M4, tau_avg_all_M4, E_val, 'linear', NaN);
        tau5 = tau_avg_all_M5_valid(j);
        
        sca1 = interp1(E_avg_all_M1, sca_avg_all_M1, E_val, 'linear', NaN);
        sca2 = interp1(E_avg_all_M2, sca_avg_all_M2, E_val, 'linear', NaN);
        sca3 = interp1(E_avg_all_M3, sca_avg_all_M3, E_val, 'linear', NaN);
        sca4 = interp1(E_avg_all_M4, sca_avg_all_M4, E_val, 'linear', NaN);
        sca5 = sca_avg_all_M5_valid(j);
        
        fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                E_val, tau1, tau2, tau3, tau4, tau5, sca1, sca2, sca3, sca4, sca5);
    end
else
    n_pts = min([length(E_avg_all_M1), length(E_avg_all_M2), length(E_avg_all_M3), length(E_avg_all_M4)]);
    for j = 1:n_pts
        fprintf(fid, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                E_avg_all_M1(j), tau_avg_all_M1(j), tau_avg_all_M2(j), ...
                tau_avg_all_M3(j), tau_avg_all_M4(j), tau_avg_all_M4(j), ...
                sca_avg_all_M1(j), sca_avg_all_M2(j), sca_avg_all_M3(j), sca_avg_all_M4(j), sca_avg_all_M4(j));
    end
end
fclose(fid);
fprintf('✓ Overall data saved to: %s\n', filename);

% SAVE PER-BAND DATA
band_names = {'Band1', 'Band2', 'Band3'};

for b = 1:numBands
    % File for relaxation times (tau)
    filename_tau = fullfile(outputDir, sprintf('%s_%s_tau_5methods.dat', baseName, band_names{b}));
    fid_tau = fopen(filename_tau, 'w');
    fprintf(fid_tau, '# %s - Relaxation times by method\n', band_names{b});
    fprintf(fid_tau, '# M5: tau = TDF.xx(E, EF=0, T=300K) / (V^2(E)*g(E))\n');
    fprintf(fid_tau, '# Energy  Tau_M1  Tau_M2  Tau_M3  Tau_M4  Tau_M5  V2g_sum  TDF_value\n');
    
    % File for scattering rates
    filename_sca = fullfile(outputDir, sprintf('%s_%s_sca_5methods.dat', baseName, band_names{b}));
    fid_sca = fopen(filename_sca, 'w');
    fprintf(fid_sca, '# %s - Scattering rates by method\n', band_names{b});
    fprintf(fid_sca, '# M5: S = 1/tau where tau = TDF.xx / (V^2*g)\n');
    fprintf(fid_sca, '# Energy  Sca_M1  Sca_M2  Sca_M3  Sca_M4  Sca_M5\n');
    
    % Check if band has data
    if ~isempty(E_avg_by_band_M1{b})
        n_band_pts = length(E_avg_by_band_M1{b});
        
        for j = 1:n_band_pts
            E_val = E_avg_by_band_M1{b}(j);
            
            % Get tau values
            tau1 = tau_avg_by_band_M1{b}(j);
            tau2 = tau_avg_by_band_M2{b}(j);
            tau3 = tau_avg_by_band_M3{b}(j);
            tau4 = tau_avg_by_band_M4{b}(j);
            
            % Get scattering rate values
            sca1 = sca_avg_by_band_M1{b}(j);
            sca2 = sca_avg_by_band_M2{b}(j);
            sca3 = sca_avg_by_band_M3{b}(j);
            sca4 = sca_avg_by_band_M4{b}(j);
            
            % Get M5 values if available
            if exist('E_avg_by_band_M5', 'var') && ~isempty(E_avg_by_band_M5{b})
                tau5 = tau_avg_by_band_M5{b}(j);
                sca5 = sca_avg_by_band_M5{b}(j);
            else
                tau5 = NaN;
                sca5 = NaN;
            end
            
            % Get V2g and TDF values for M5 if available
            V2g_val = NaN;
            TDF_val = NaN;
            if exist('V2g_avg_by_band', 'var') && ~isempty(V2g_avg_by_band{b})
                if j <= length(V2g_avg_by_band{b})
                    V2g_val = V2g_avg_by_band{b}(j);
                end
            end
            if exist('TDF_avg_by_band', 'var') && ~isempty(TDF_avg_by_band{b})
                if j <= length(TDF_avg_by_band{b})
                    TDF_val = TDF_avg_by_band{b}(j);
                end
            end
            
            % Write tau data
            fprintf(fid_tau, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                    E_val, tau1, tau2, tau3, tau4, tau5, V2g_val, TDF_val);
            
            % Write scattering rate data
            fprintf(fid_sca, '%.6f  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                    E_val, sca1, sca2, sca3, sca4, sca5);
        end
    end
    
    fclose(fid_tau);
    fclose(fid_sca);
    fprintf('✓ %s tau data saved to: %s\n', band_names{b}, filename_tau);
    fprintf('✓ %s scattering data saved to: %s\n', band_names{b}, filename_sca);
end

% SAVE INDIVIDUAL K-POINT DATA PER BAND
fprintf('\n--- Saving individual k-point data ---\n');
for b = 1:numBands
    if ~isempty(all_E_by_band{b})
        filename_kpt = fullfile(outputDir, sprintf('%s_%s_kpoints.dat', baseName, band_names{b}));
        fid_kpt = fopen(filename_kpt, 'w');
        fprintf(fid_kpt, '# %s - Individual k-point data\n', band_names{b});
        fprintf(fid_kpt, '# Energy  Scattering_Rate  Tau  DOS  V^2\n');
        
        n_kpts = length(all_E_by_band{b});
        for k = 1:n_kpts
            fprintf(fid_kpt, '%.6f  %.6e  %.6e  %.6e  %.6e\n', ...
                    all_E_by_band{b}(k), all_sca_by_band{b}(k), ...
                    all_tau_by_band{b}(k), all_DOS_by_band{b}(k), ...
                    all_V2_by_band{b}(k));
        end
        
        fclose(fid_kpt);
        fprintf('✓ %s k-points data saved to: %s (%d points)\n', ...
                band_names{b}, filename_kpt, n_kpts);
    end
end

% SAVE COMBINED K-POINT DATA (all bands)
if ~isempty(all_E)
    filename_all_kpt = fullfile(outputDir, sprintf('%s_all_kpoints.dat', baseName));
    fid_all_kpt = fopen(filename_all_kpt, 'w');
    fprintf(fid_all_kpt, '# All bands combined - Individual k-point data\n');
    fprintf(fid_all_kpt, '# Band  Energy  Scattering_Rate  Tau  DOS  V^2\n');
    
    for k = 1:length(all_E)
        fprintf(fid_all_kpt, '%d  %.6f  %.6e  %.6e  %.6e  %.6e\n', ...
                all_band(k), all_E(k), all_sca(k), all_tau(k), ...
                all_DOS(k), all_V2(k));
    end
    
    fclose(fid_all_kpt);
    fprintf('✓ Combined k-points data saved to: %s (%d points)\n', ...
            filename_all_kpt, length(all_E));
end

% SAVE SUMMARY STATISTICS FILE
filename_stats = fullfile(outputDir, sprintf('%s_summary_statistics.dat', baseName));
fid_stats = fopen(filename_stats, 'w');
fprintf(fid_stats, '# Summary Statistics for All Methods\n');
fprintf(fid_stats, '# \n');
fprintf(fid_stats, '# Overall Statistics:\n');
fprintf(fid_stats, '# Method  Mean_Tau  Std_Tau  Min_Tau  Max_Tau  Mean_Sca\n');

methods_tau = {tau_avg_all_M1, tau_avg_all_M2, tau_avg_all_M3, tau_avg_all_M4};
methods_sca = {sca_avg_all_M1, sca_avg_all_M2, sca_avg_all_M3, sca_avg_all_M4};
methods_names = {'M1: Simple', 'M2: DOS-wt S', 'M3: DOS-wt tau', 'M4: V2g-wt tau'};

for m = 1:4
    valid_tau = methods_tau{m}(~isnan(methods_tau{m}) & methods_tau{m} > 0);
    valid_sca = methods_sca{m}(~isnan(methods_sca{m}) & methods_sca{m} > 0);
    if ~isempty(valid_tau)
        fprintf(fid_stats, '%s  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
                methods_names{m}, mean(valid_tau), std(valid_tau), ...
                min(valid_tau), max(valid_tau), mean(valid_sca));
    end
end

% M5 statistics
if exist('tau_avg_all_M5_valid', 'var')
    fprintf(fid_stats, 'M5: TDF/V2g  %.6e  %.6e  %.6e  %.6e  %.6e\n', ...
            mean(tau_avg_all_M5_valid), std(tau_avg_all_M5_valid), ...
            min(tau_avg_all_M5_valid), max(tau_avg_all_M5_valid), ...
            mean(sca_avg_all_M5_valid));
end

% Per-band statistics
for b = 1:numBands
    fprintf(fid_stats, '\n# Band %d Statistics:\n', b);
    fprintf(fid_stats, '# Method  Mean_Tau  Std_Tau  Mean_Sca\n');
    
    methods_tau_band = {tau_avg_by_band_M1{b}, tau_avg_by_band_M2{b}, ...
                        tau_avg_by_band_M3{b}, tau_avg_by_band_M4{b}};
    methods_sca_band = {sca_avg_by_band_M1{b}, sca_avg_by_band_M2{b}, ...
                        sca_avg_by_band_M3{b}, sca_avg_by_band_M4{b}};
    
    for m = 1:4
        if ~isempty(methods_tau_band{m})
            valid_tau = methods_tau_band{m}(~isnan(methods_tau_band{m}) & methods_tau_band{m} > 0);
            valid_sca = methods_sca_band{m}(~isnan(methods_sca_band{m}) & methods_sca_band{m} > 0);
            if ~isempty(valid_tau)
                fprintf(fid_stats, '%s  %.6e  %.6e  %.6e\n', ...
                        methods_names{m}, mean(valid_tau), std(valid_tau), mean(valid_sca));
            end
        end
    end
    
    % M5 per-band
    if exist('tau_avg_by_band_M5', 'var') && ~isempty(tau_avg_by_band_M5{b})
        valid_tau5 = tau_avg_by_band_M5{b}(~isnan(tau_avg_by_band_M5{b}) & tau_avg_by_band_M5{b} > 0);
        valid_sca5 = sca_avg_by_band_M5{b}(~isnan(sca_avg_by_band_M5{b}) & sca_avg_by_band_M5{b} > 0);
        if ~isempty(valid_tau5)
            fprintf(fid_stats, 'M5: TDF/V2g  %.6e  %.6e  %.6e\n', ...
                    mean(valid_tau5), std(valid_tau5), mean(valid_sca5));
        end
    end
end

fclose(fid_stats);
fprintf('✓ Summary statistics saved to: %s\n', filename_stats);

fprintf('\n=== All data saved successfully ===\n');
fprintf('Output directory: %s\n', outputDir);
fprintf('Files created:\n');
fprintf('  - *_overall_5methods.dat (combined bands)\n');
fprintf('  - *_Band1_2_3_tau_5methods.dat (per-band tau)\n');
fprintf('  - *_Band1_2_3_sca_5methods.dat (per-band scattering)\n');
fprintf('  - *_Band1_2_3_kpoints.dat (individual k-points)\n');
fprintf('  - *_all_kpoints.dat (all k-points combined)\n');
fprintf('  - *_summary_statistics.dat (statistics)\n');
fprintf('\n=== Complete ===\n');