%% First, let's check the TDF_n and TDF_ph_n structures
fprintf('Checking band-dependent TDF variables...\n\n');

% Check TDF_n (total band-resolved TDF)
if exist('TDF_n', 'var')
    fprintf('✓ Found TDF_n\n');
    if isfield(TDF_n, 'xx')
        fprintf('  TDF_n.xx size: [%d x %d x %d x %d]\n', size(TDF_n.xx));
        fprintf('  Dim 1: Energy (%d points)\n', size(TDF_n.xx, 1));
        fprintf('  Dim 2: Bands (%d bands)\n', size(TDF_n.xx, 2));
        fprintf('  Dim 3: Fermi levels (%d)\n', size(TDF_n.xx, 3));
        fprintf('  Dim 4: Temperatures (%d)\n', size(TDF_n.xx, 4));
    end
    % Check other components
    fields_n = fieldnames(TDF_n);
    fprintf('  Available components: %s\n', strjoin(fields_n, ', '));
else
    fprintf('✗ TDF_n not found\n');
end

% Check TDF_ph_n (phonon band-resolved TDF)
if exist('TDF_ph_n', 'var')
    fprintf('\n✓ Found TDF_ph_n\n');
    if isfield(TDF_ph_n, 'xx')
        fprintf('  TDF_ph_n.xx size: [%d x %d x %d x %d]\n', size(TDF_ph_n.xx));
    end
    fields_ph = fieldnames(TDF_ph_n);
    fprintf('  Available components: %s\n', strjoin(fields_ph, ', '));
else
    fprintf('\n✗ TDF_ph_n not found\n');
end

% Check TDF (total TDF, summed over bands)
if exist('TDF', 'var')
    fprintf('\n✓ Found TDF (total)\n');
    if isfield(TDF, 'xx')
        fprintf('  TDF.xx size: [%d x %d x %d]\n', size(TDF.xx));
    end
end