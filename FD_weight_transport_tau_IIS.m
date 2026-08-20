%% Load the data files
load 'TE_ZrNiSn_kScan_holes.mat';
load 'RelaxTimes_IIS_ZrNiSn_kScan_holes.mat';

%% Initialize
figure1 = figure('Name', 'Scattering Rates');
fig1 = gcf; fig1.Position(3:4) = [1000, 800]; hold on;

figure2 = figure('Name', 'Relaxation Times');
fig2 = gcf; fig2.Position(3:4) = [1000, 800]; hold on;

numBands = 3;
palette = [
    0.1216 0.4667 0.7059;
    1.0000 0.4980 0.0549;
    0.1725 0.6275 0.1725;
    0.8392 0.1529 0.1569;
    0.5804 0.4039 0.7412;
    0.5490 0.3373 0.2941;
    0.8902 0.4667 0.7608;
    0.4980 0.4980 0.4980;
];

% Initialize overall arrays
E_avg_all_M1=[]; sca_avg_all_M1=[]; tau_avg_all_M1=[];
E_avg_all_M2=[]; sca_avg_all_M2=[]; tau_avg_all_M2=[];
E_avg_all_M3=[]; sca_avg_all_M3=[]; tau_avg_all_M3=[];
E_avg_all_M4=[]; sca_avg_all_M4=[]; tau_avg_all_M4=[];
E_avg_all_M5=[]; sca_avg_all_M5=[]; tau_avg_all_M5=[];
E_avg_all_V2g=[]; V2g_all_total=[];

% Initialize per-band cell arrays
[E_avg_by_band_M1, sca_avg_by_band_M1, tau_avg_by_band_M1] = deal(cell(numBands,1));
[E_avg_by_band_M2, sca_avg_by_band_M2, tau_avg_by_band_M2] = deal(cell(numBands,1));
[E_avg_by_band_M3, sca_avg_by_band_M3, tau_avg_by_band_M3] = deal(cell(numBands,1));
[E_avg_by_band_M4, sca_avg_by_band_M4, tau_avg_by_band_M4] = deal(cell(numBands,1));
[E_avg_by_band_M5, sca_avg_by_band_M5, tau_avg_by_band_M5] = deal(cell(numBands,1));
[E_avg_by_band_V2g, V2g_avg_by_band] = deal(cell(numBands,1));

% Individual k-points
[all_E_by_band, all_sca_by_band, all_tau_by_band, all_DOS_by_band, all_V2_by_band] = deal(cell(numBands,1));
all_E=[]; all_sca=[]; all_tau=[]; all_band=[]; all_DOS=[]; all_V2=[];

% Configuration
EF_index=8; T_index=1; T_temperature=300;
k_B=8.617333262145e-5; E_F=0;

% Fermi window function
fermi_window = @(E) (1/(k_B*T_temperature)) * exp((E-E_F)/(k_B*T_temperature)) ./ ...
                     (1 + exp((E-E_F)/(k_B*T_temperature))).^2;

% Check TDF
if exist('TDF','var') && isfield(TDF,'xx')
    fprintf('TDF.xx found. Dims: [%d,%d,%d]\n', size(TDF.xx));
    use_TDF=true;
else
    fprintf('TDF not found.\n');
    use_TDF=false;
end

% Energy range
if exist('E_array','var')
    E_min=min(E_array); E_max=max(E_array);
    fprintf('E_array: %.4f to %.4f eV (%d points)\n', E_min, E_max, length(E_array));
    energy_indices=find(arrayfun(@(i) ~isempty(state_ID(i,1).E) && state_ID(i,1).E>=E_min && state_ID(i,1).E<=E_max, 1:size(state_ID,1)));
    max_idx=length(energy_indices);
else
    max_idx=size(state_ID,1); energy_indices=1:max_idx;
    E_min=0; E_max=0.5;
end

fprintf('Processing %d energy points...\n', max_idx);

%% Main processing loop
for idx = 1:max_idx
    i = energy_indices(idx);
    E_row=[]; sca_M1_row=[]; tau_M1_row=[]; sca_M2_row=[]; tau_M2_row=[];
    sca_M3_row=[]; tau_M3_row=[]; sca_M4_row=[]; tau_M4_row=[];
    DOS_row=[]; V2g_band_row=[];
    
    for b = 1:numBands
        if ~isempty(taus_IIS(i,b).x) && size(taus_IIS(i,b).x,3)>=1
            E = state_ID(i,b).E;
            tau_3D = taus_IIS(i,b).x;
            
            if size(tau_3D,1)>=EF_index && size(tau_3D,2)>=T_index
                tau_k = squeeze(tau_3D(EF_index, T_index, :));
                
                V2_k = ones(size(tau_k));
                if isfield(state_ID(i,b),'V')
                    V_k = state_ID(i,b).V(:);
                    if length(V_k)==length(tau_k), V2_k=V_k.^2; end
                end
                
                DOS_k = ones(size(tau_k));
                if isfield(state_ID(i,b),'DOS') && ~isempty(state_ID(i,b).DOS)
                    dk = state_ID(i,b).DOS(:);
                    if length(dk)==length(tau_k), DOS_k=dk; end
                end
                
                pos = tau_k>0 & V2_k>0;
                
                if any(pos)
                    tk = tau_k(pos); v2 = V2_k(pos); dk = DOS_k(pos);
                    sk = 1./tk; npts = length(sk);
                    
                    % Store individual k-points
                    all_E_by_band{b}=[all_E_by_band{b}; repmat(E,npts,1)];
                    all_sca_by_band{b}=[all_sca_by_band{b}; sk];
                    all_tau_by_band{b}=[all_tau_by_band{b}; tk];
                    all_DOS_by_band{b}=[all_DOS_by_band{b}; dk];
                    all_V2_by_band{b}=[all_V2_by_band{b}; v2];
                    
                    all_E=[all_E; repmat(E,npts,1)];
                    all_sca=[all_sca; sk]; all_tau=[all_tau; tk];
                    all_band=[all_band; repmat(b,npts,1)];
                    all_DOS=[all_DOS; dk]; all_V2=[all_V2; v2];
                    
                    % M1: Simple
                    s1=mean(sk); t1=1/s1;
                    
                    % M2: DOS-weighted S
                    s2 = sum(sk.*dk)/sum(dk); t2=1/s2;
                    
                    % M3: DOS-weighted tau
                    t3 = sum(tk.*dk)/sum(dk); s3=1/t3;
                    
                    % M4: Transport (V^2*g)
                    w4 = v2.*dk;
                    t4 = sum(tk.*w4)/sum(w4); s4=1/t4;
                    
                    v2g = sum(v2.*dk);
                    
                    % Store per-band
                    E_avg_by_band_M1{b}(end+1)=E; sca_avg_by_band_M1{b}(end+1)=s1; tau_avg_by_band_M1{b}(end+1)=t1;
                    E_avg_by_band_M2{b}(end+1)=E; sca_avg_by_band_M2{b}(end+1)=s2; tau_avg_by_band_M2{b}(end+1)=t2;
                    E_avg_by_band_M3{b}(end+1)=E; sca_avg_by_band_M3{b}(end+1)=s3; tau_avg_by_band_M3{b}(end+1)=t3;
                    E_avg_by_band_M4{b}(end+1)=E; sca_avg_by_band_M4{b}(end+1)=s4; tau_avg_by_band_M4{b}(end+1)=t4;
                    E_avg_by_band_V2g{b}(end+1)=E; V2g_avg_by_band{b}(end+1)=v2g;
                    
                    % Overall rows
                    E_row(end+1)=E;
                    sca_M1_row(end+1)=s1; tau_M1_row(end+1)=t1;
                    sca_M2_row(end+1)=s2; tau_M2_row(end+1)=t2;
                    sca_M3_row(end+1)=s3; tau_M3_row(end+1)=t3;
                    sca_M4_row(end+1)=s4; tau_M4_row(end+1)=t4;
                    DOS_row(end+1)=sum(dk);
                    V2g_band_row(end+1)=v2g;
                end
            end
        end
    end
    
    if ~isempty(E_row) && sum(DOS_row)>0
        E_avg = sum(E_row.*DOS_row)/sum(DOS_row);
        V2g_total = sum(V2g_band_row);
        
        E_avg_all_M1(end+1)=E_avg; sca_avg_all_M1(end+1)=sum(sca_M1_row.*DOS_row)/sum(DOS_row); tau_avg_all_M1(end+1)=sum(tau_M1_row.*DOS_row)/sum(DOS_row);
        E_avg_all_M2(end+1)=E_avg; sca_avg_all_M2(end+1)=sum(sca_M2_row.*DOS_row)/sum(DOS_row); tau_avg_all_M2(end+1)=sum(tau_M2_row.*DOS_row)/sum(DOS_row);
        E_avg_all_M3(end+1)=E_avg; sca_avg_all_M3(end+1)=sum(sca_M3_row.*DOS_row)/sum(DOS_row); tau_avg_all_M3(end+1)=sum(tau_M3_row.*DOS_row)/sum(DOS_row);
        E_avg_all_M4(end+1)=E_avg; sca_avg_all_M4(end+1)=sum(sca_M4_row.*DOS_row)/sum(DOS_row); tau_avg_all_M4(end+1)=sum(tau_M4_row.*DOS_row)/sum(DOS_row);
        E_avg_all_V2g(end+1)=E_avg; V2g_all_total(end+1)=V2g_total;
    end
    
    if mod(idx,50)==0, fprintf('  %d/%d\n', idx, max_idx); end
end

N = length(E_avg_all_M1);
fprintf('Done. %d overall energy points.\n', N);

%% M5: TDF/(V^2*g)
fprintf('Calculating M5...\n');

% Initialize per-band M5
for b=1:numBands
    nb=length(E_avg_by_band_V2g{b});
    E_avg_by_band_M5{b}=zeros(nb,1);
    tau_avg_by_band_M5{b}=zeros(nb,1);
    sca_avg_by_band_M5{b}=zeros(nb,1);
end

E_avg_all_M5=zeros(N,1); tau_avg_all_M5=zeros(N,1); sca_avg_all_M5=zeros(N,1);

if use_TDF
    TDF_xx = TDF.xx(:, EF_index, T_index);
    
    for j=1:N
        E_val=E_avg_all_M1(j);
        V2g_val=V2g_all_total(j);
        TDF_val=interp1(E_array, TDF_xx, E_val, 'linear', 0);
        E_avg_all_M5(j)=E_val;
        if V2g_val>0 && TDF_val>0
            tau_avg_all_M5(j)=TDF_val/V2g_val;
            sca_avg_all_M5(j)=1/tau_avg_all_M5(j);
        else
            tau_avg_all_M5(j)=NaN; sca_avg_all_M5(j)=NaN;
        end
    end
    
    for b=1:numBands
        nb=length(E_avg_by_band_V2g{b});
        for j=1:nb
            E_val=E_avg_by_band_V2g{b}(j);
            V2g_val=V2g_avg_by_band{b}(j);
            TDF_val=interp1(E_array, TDF_xx, E_val, 'linear', 0);
            E_avg_by_band_M5{b}(j)=E_val;
            if V2g_val>0 && TDF_val>0
                tau_avg_by_band_M5{b}(j)=TDF_val/V2g_val;
                sca_avg_by_band_M5{b}(j)=1/tau_avg_by_band_M5{b}(j);
            else
                tau_avg_by_band_M5{b}(j)=NaN; sca_avg_by_band_M5{b}(j)=NaN;
            end
        end
    end
else
    E_avg_all_M5=E_avg_all_M4; tau_avg_all_M5=tau_avg_all_M4; sca_avg_all_M5=sca_avg_all_M4;
    for b=1:numBands
        E_avg_by_band_M5{b}=E_avg_by_band_M4{b};
        tau_avg_by_band_M5{b}=tau_avg_by_band_M4{b};
        sca_avg_by_band_M5{b}=sca_avg_by_band_M4{b};
    end
end

%% FD Integration - Overall
fprintf('\n=== FD Integration (T=%dK) ===\n', T_temperature);

E_vec = E_avg_all_M1(:);
fw = fermi_window(E_vec);
V2g_vec = V2g_all_total(:);
weight = V2g_vec .* fw;
wsum = sum(weight);

tau_FD_M1 = sum(tau_avg_all_M1(:).*weight)/wsum;
tau_FD_M2 = sum(tau_avg_all_M2(:).*weight)/wsum;
tau_FD_M3 = sum(tau_avg_all_M3(:).*weight)/wsum;
tau_FD_M4 = sum(tau_avg_all_M4(:).*weight)/wsum;

v5 = ~isnan(tau_avg_all_M5) & tau_avg_all_M5>0;
if sum(v5)>0
    tau_FD_M5 = sum(tau_avg_all_M5(v5).*weight(v5))/sum(weight(v5));
else
    tau_FD_M5=NaN;
end

fprintf('\nOverall FD-Integrated Tau:\n');
fprintf('  M1: %.4f fs\n', tau_FD_M1);
fprintf('  M2: %.4f fs\n', tau_FD_M2);
fprintf('  M3: %.4f fs\n', tau_FD_M3);
fprintf('  M4: %.4f fs\n', tau_FD_M4);
fprintf('  M5: %.4f fs\n', tau_FD_M5);

%% FD Integration - Per Band
fprintf('\n=== Per-Band FD Integration ===\n');

tau_FD_band = zeros(5, numBands);
tau_avg_band = zeros(5, numBands);

for b = 1:numBands
    if isempty(E_avg_by_band_M1{b}), continue; end
    
    % Get per-band data as column vectors
    Eb = E_avg_by_band_M1{b}(:);
    V2gb = V2g_avg_by_band{b}(:);
    fwb = fermi_window(Eb);
    
    % Make sure all have same length
    n_pts = length(Eb);
    
    % Get tau arrays
    t1 = tau_avg_by_band_M1{b}(:);
    t2 = tau_avg_by_band_M2{b}(:);
    t3 = tau_avg_by_band_M3{b}(:);
    t4 = tau_avg_by_band_M4{b}(:);
    t5 = tau_avg_by_band_M5{b}(:);
    
    % Create valid mask (same length for all)
    valid = true(n_pts,1);
    valid = valid & ~isnan(t1) & t1>0;
    valid = valid & V2gb>0 & fwb>0;
    
    if sum(valid) == 0, continue; end
    
    wb = V2gb(valid) .* fwb(valid);
    wbs = sum(wb);
    
    % Calculate FD tau for each method
    tau_list = {t1, t2, t3, t4, t5};
    for m = 1:5
        tm = tau_list{m};
        if length(tm) == n_pts
            tm_valid = tm(valid);
            vm = ~isnan(tm_valid) & tm_valid > 0;
            if sum(vm) > 0
                tau_FD_band(m,b) = sum(tm_valid(vm).*wb(vm))/sum(wb(vm));
                tau_avg_band(m,b) = mean(tm_valid(vm));
            end
        end
    end
    
    fprintf('\nBand %d:\n', b);
    fprintf('  Method          FD-Tau[fs]    Avg-Tau[fs]\n');
    names = {'M1-Simple','M2-DOS-S','M3-DOS-tau','M4-V2g','M5-TDF'};
    for m = 1:5
        if tau_FD_band(m,b) > 0
            fprintf('  %-12s    %8.4f      %8.4f\n', names{m}, tau_FD_band(m,b), tau_avg_band(m,b));
        end
    end
end

%% Summary Table
fprintf('\n========================================\n');
fprintf('FINAL SUMMARY TABLE [fs]\n');
fprintf('========================================\n');
fprintf('%-12s %8s %8s %8s %8s\n', 'Method', 'Band1', 'Band2', 'Band3', 'Overall');
names = {'M1-Simple','M2-DOS-S','M3-DOS-tau','M4-V2g','M5-TDF'};
overall = [tau_FD_M1, tau_FD_M2, tau_FD_M3, tau_FD_M4, tau_FD_M5];
for m = 1:5
    fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', names{m}, ...
            tau_FD_band(m,1), tau_FD_band(m,2), tau_FD_band(m,3), overall(m));
end

%% Figures
% Figure 1: Scattering Rates
figure(figure1);
valid5 = ~isnan(tau_avg_all_M5) & tau_avg_all_M5>0;
if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_sca_by_band{1}, 5, [0.7 0.7 0.7], 'filled'); hold on;
end
plot(E_avg_all_M1, sca_avg_all_M1, 'k-', 'LineWidth', 2);
plot(E_avg_all_M2, sca_avg_all_M2, 'b--', 'LineWidth', 2);
plot(E_avg_all_M3, sca_avg_all_M3, 'g-.', 'LineWidth', 2);
plot(E_avg_all_M4, sca_avg_all_M4, 'r:', 'LineWidth', 2);
if sum(valid5)>0
    plot(E_avg_all_M5(valid5), sca_avg_all_M5(valid5), 'm-', 'LineWidth', 2.5);
end
set(gca, 'YScale', 'log'); xlim([E_min, E_max]); grid on;
xlabel('Energy [eV]'); ylabel('Scattering Rate [fs^{-1}]');
title('Scattering Rates'); legend({'k-pts','M1','M2','M3','M4','M5'}, 'Location','best');

% Figure 2: Relaxation Times
figure(figure2);
yyaxis right;
fw_norm = fw/max(fw);
area(E_vec, fw_norm, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
ylabel('Fermi Window'); set(gca,'YColor',[0.5 0.5 0.5]); ylim([0,1.1]);

yyaxis left;
if ~isempty(all_E_by_band{1})
    scatter(all_E_by_band{1}, all_tau_by_band{1}, 5, [0.7 0.7 0.7], 'filled'); hold on;
end
plot(E_avg_all_M1, tau_avg_all_M1, 'k-', 'LineWidth', 2);
plot(E_avg_all_M2, tau_avg_all_M2, 'b--', 'LineWidth', 2);
plot(E_avg_all_M3, tau_avg_all_M3, 'g-.', 'LineWidth', 2);
plot(E_avg_all_M4, tau_avg_all_M4, 'r:', 'LineWidth', 2);
if sum(valid5)>0
    plot(E_avg_all_M5(valid5), tau_avg_all_M5(valid5), 'm-', 'LineWidth', 2.5);
end
set(gca, 'YScale', 'log'); xlim([E_min, E_max]); grid on;
xlabel('Energy [eV]'); ylabel('\tau [fs]');
title(sprintf('Relaxation Times (T=%dK)', T_temperature));
legend({'k-pts','M1','M2','M3','M4','M5','Fermi'}, 'Location','best');

% Figure 3: Bar chart
figure('Name','FD Results');
tau_vals = [tau_FD_M1, tau_FD_M2, tau_FD_M3, tau_FD_M4, tau_FD_M5];
bar(tau_vals);
set(gca, 'XTickLabel', {'M1','M2','M3','M4','M5'});
ylabel('\tau_{FD} [fs]'); title(sprintf('FD-Integrated Tau (T=%dK)', T_temperature));
grid on;
for k=1:5
    if ~isnan(tau_vals(k))
        text(k, tau_vals(k)+0.02*max(tau_vals), sprintf('%.4f',tau_vals(k)), ...
             'HorizontalAlignment','center');
    end
end

fprintf('\n=== Complete ===\n');