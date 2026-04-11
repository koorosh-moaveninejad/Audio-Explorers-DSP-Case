clear; clc; close all;

% 1. Configuration & Paths ---
basePath = '/Users/kourosh/Desktop/University/Self Study/Audio-Explorers20206/Diagnostics DSP/Patient Data';
templatePath = fullfile(basePath, 'lostOaes.json');

% Load Templates and Folders
templates = jsondecode(fileread(templatePath));
temp_names = fieldnames(templates);

dirInfo = dir(fullfile(basePath, 'patient_*'));
patientFolders = {dirInfo.name};

numPatients = length(patientFolders);

% Pre-allocate storage
all_scores_matrix = []; 
oae_results_cell = cell(numPatients, 1); 
snr_values = zeros(numPatients, 1);
repeat_corr_values = zeros(numPatients, 1);
fs_values = zeros(numPatients, 1);

fprintf('Mission Start: Recovering %d newborns...\n', numPatients);

% --- Optional: Template self-correlation analysis ---
numTemps = length(temp_names);
template_corr = zeros(numTemps, numTemps);

for i = 1:numTemps
    t1 = templates.(temp_names{i})(:);
    t1 = t1 / (norm(t1) + eps);

    for j = 1:numTemps
        t2 = templates.(temp_names{j})(:);
        t2_rs = resample(t2, length(t1), length(t2));
        t2_rs = t2_rs / (norm(t2_rs) + eps);
        template_corr(i,j) = dot(t1, t2_rs);
    end
end

disp('--- TEMPLATE SELF-CORRELATION MATRIX ---');
disp(array2table(template_corr, 'VariableNames', temp_names, 'RowNames', temp_names));

figure('Name','Template Self-Correlation','Color','w');
imagesc(template_corr);
colorbar;
axis square;
title('Template Self-Correlation Matrix');
xticks(1:numTemps); xticklabels(temp_names); xtickangle(45);
yticks(1:numTemps); yticklabels(temp_names);

fprintf('\n--- MAX TEMPLATE CONFUSION ---\n');
for i = 1:numTemps
    vals = template_corr(i,:);
    vals(i) = -inf;
    fprintf('%s -> max confusion = %.3f\n', temp_names{i}, max(vals));
end

% --- 2. Main Processing Loop ---
for p = 1:numPatients
    pID = strrep(patientFolders{p}, 'patient_', '');
    pFolder = fullfile(basePath, patientFolders{p});
    
    % A. Load Data & Metadata
    [rec, fs] = audioread(fullfile(pFolder, ['Patient_' pID '_rec.wav']));
    fs_values(p) = fs;
    info = jsondecode(fileread(fullfile(pFolder, ['Patient_' pID '_info.json'])));
    
    % Read target frequencies robustly
    if isfield(info, 'centerFreq')
        target_freqs = info.centerFreq(:);
    elseif isfield(info, 'centerFreqs')
        target_freqs = info.centerFreqs(:);
    else
        target_freqs = [1000; 1414.2; 2000; 2828.4; 4000];
    end
    
    % B. Extraction & Sub-Sample Alignment
    types = {'A','B','C','D'};
    epochData = cell(1,4);
    p_ref = [];
    
    for i = 1:4
        act = audioread(fullfile(pFolder, ['Patient_' pID '_active' types{i} '.wav']));
        starts = find(diff([0; act > 0.5]) == 1);
        tmp = [];
        
        for s = starts'
            e = s + info.epochSize - 1;
            if e <= length(rec)
                seg = rec(s:e);
                seg = seg - mean(seg);
                
                % Jitter correction via cross-correlation on stimulus onset
                if isempty(p_ref)
                    p_ref = seg(1:min(100, length(seg)));
                    shifted = seg;
                else
                    seg_head = seg(1:min(100, length(seg)));
                    ref_head = p_ref(1:min(length(p_ref), length(seg_head)));
                    
                    [c_a, lags_a] = xcorr(seg_head(1:length(ref_head)), ref_head, 15, 'coeff');
                    [~, m_idx_align] = max(c_a);
                    shifted = circshift(seg, -lags_a(m_idx_align));
                end
                
                tmp = [tmp, shifted]; %#ok<AGROW>
            end
        end
        
        epochData{i} = tmp;
    end
    
    % Safety check
    minE = min(cellfun(@(x) size(x,2), epochData));
    if minE == 0
        warning('Patient %s has no valid epochs. Marking as REFER-like by default.', pID);
        oae_results_cell{p} = zeros(info.epochSize, 1);
        snr_values(p) = -Inf;
        repeat_corr_values(p) = 0;
        
        for j = 1:length(temp_names)
            all_scores_matrix = [all_scores_matrix; p, j, 0]; %#ok<AGROW>
        end
        continue;
    end
    
    % C. Nonlinear Summation (+1+1+1-3)
    acc = zeros(info.epochSize, 1);
    acc1 = zeros(info.epochSize, 1);
    acc2 = zeros(info.epochSize, 1);
    v_sets = 0;
    v1 = 0;
    v2 = 0;
    
    for k = 1:minE
        quad = epochData{1}(:,k) + epochData{2}(:,k) + epochData{3}(:,k) + epochData{4}(:,k);
        
        if all(isfinite(quad))
            acc = acc + quad;
            v_sets = v_sets + 1;
            
            if mod(k,2) == 1
                acc1 = acc1 + quad;
                v1 = v1 + 1;
            else
                acc2 = acc2 + quad;
                v2 = v2 + 1;
            end
        end
    end
    
    oae_raw  = acc  / max(1, v_sets);
    oae_raw1 = acc1 / max(1, v1);
    oae_raw2 = acc2 / max(1, v2);
    
    % D. Post-Processing (4-16 ms gating + bandpass)
    t = (0:info.epochSize-1)' / fs;
    win_env = exp(-((t - 0.008).^2) / (2 * 0.004^2)); 
    
    bpFilt = designfilt('bandpassiir', 'FilterOrder', 10, ...
        'HalfPowerFrequency1', 1000, ...
        'HalfPowerFrequency2', 3800, ...
        'SampleRate', fs);
    
    oae_clean  = filtfilt(bpFilt, oae_raw  .* win_env);
    oae_clean1 = filtfilt(bpFilt, oae_raw1 .* win_env);
    oae_clean2 = filtfilt(bpFilt, oae_raw2 .* win_env);
    
    % E. Spectral Denoising / noise subtraction
    noise_est = (mean(epochData{1},2) + mean(epochData{2},2)) - ...
                (mean(epochData{3},2) + mean(epochData{4},2));
    noise_filt = filtfilt(bpFilt, noise_est);
    
    oae_clean  = oae_clean  - 0.1 * noise_filt;
    oae_clean1 = oae_clean1 - 0.1 * noise_filt;
    oae_clean2 = oae_clean2 - 0.1 * noise_filt;
    
    oae_results_cell{p} = oae_clean;
    
    % F. SNR & Repeatability Metrics
    resp_idx  = (t > 0.004 & t < 0.016);
    noise_idx = (t > 0.017 & t < 0.020);
    
    resp_rms = rms(oae_clean(resp_idx));
    noise_rms = rms(oae_clean(noise_idx));
    snr_val = 20 * log10(resp_rms / max(noise_rms, eps));
    snr_values(p) = snr_val;
    
    x1 = oae_clean1(resp_idx);
    x2 = oae_clean2(resp_idx);
    
    if norm(x1) > 0 && norm(x2) > 0
        r_corr = corr(x1, x2);
    else
        r_corr = 0;
    end
    if isnan(r_corr)
        r_corr = 0;
    end
    repeat_corr_values(p) = max(0, r_corr);
    
    % G. Dual-Domain Matching
    oae_crop = oae_clean(resp_idx);
    Nfft = length(oae_crop);
    
    if Nfft < 2
        for j = 1:length(temp_names)
            all_scores_matrix = [all_scores_matrix; p, j, 0]; %#ok<AGROW>
        end
        continue;
    end
    
    Ye = abs(fft(oae_crop) / Nfft);
    Ye_half = Ye(1:floor(Nfft/2)+1);
    Ye_norm = Ye_half / (max(Ye_half) + eps);
    f_axis = fs * (0:floor(Nfft/2)) / Nfft;
    
    for j = 1:length(temp_names)
        target = templates.(temp_names{j})(:);
        
        target_adj = [target; zeros(max(0, length(oae_clean) - length(target)), 1)];
        target_adj = target_adj(1:length(oae_clean));
        target_crop = target_adj(resp_idx);
        
        % Time-domain score
        [c, ~] = xcorr(oae_crop, target_crop, 'coeff');
        time_score = max(c);
        if isnan(time_score)
            time_score = 0;
        end
        
        % Spectral score at diagnostic frequencies
        spectral_score = 0;
        Yt = abs(fft(target_crop) / Nfft);
        Yt_half = Yt(1:floor(Nfft/2)+1);
        Yt_norm = Yt_half / (max(Yt_half) + eps);
        
        for f_target = target_freqs'
            [~, bin] = min(abs(f_axis - f_target));
            spectral_score = spectral_score + (Ye_norm(bin) * Yt_norm(bin));
        end
        spectral_score = spectral_score / length(target_freqs);
        
        % Final combined similarity
        combined = (0.7 * max(0, time_score)) + (0.3 * spectral_score);
        
        all_scores_matrix = [all_scores_matrix; p, j, combined]; %#ok<AGROW>
    end
end

% --- 3. Final Result Generation ---
[~, rep_order] = sort(repeat_corr_values, 'descend');
numPass = min(8, numPatients);
pass_indices = rep_order(1:numPass);

final_mapping = table('Size', [numPatients 5], ...
    'VariableTypes', {'string', 'string', 'string', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'Result', 'Template', 'Confidence', 'RepeatCorr'});

for i = 1:numPatients
    final_mapping.PatientID(i) = string(strrep(patientFolders{i}, 'patient_', ''));
    final_mapping.Result(i) = "REFER";
    final_mapping.Template(i) = "N/A";
    final_mapping.Confidence(i) = 0;
    final_mapping.RepeatCorr(i) = repeat_corr_values(i);
end

used_templates = [];

for p_idx = pass_indices'
    p_rows = all_scores_matrix(all_scores_matrix(:,1) == p_idx, :);
    [~, s_idx] = sort(p_rows(:,3), 'descend');
    
    assigned = false;
    for r = s_idx'
        t_idx = p_rows(r,2);
        if ~any(used_templates == t_idx)
            row = find(final_mapping.PatientID == string(strrep(patientFolders{p_idx}, 'patient_', '')));
            final_mapping.Result(row) = "PASS";
            final_mapping.Template(row) = string(temp_names{t_idx});
            final_mapping.Confidence(row) = p_rows(r,3) * 100;
            used_templates = [used_templates; t_idx]; %#ok<AGROW>
            assigned = true;
            break;
        end
    end
    
    if ~assigned
        row = find(final_mapping.PatientID == string(strrep(patientFolders{p_idx}, 'patient_', '')));
        final_mapping.Result(row) = "PASS";
        final_mapping.Template(row) = "N/A";
        final_mapping.Confidence(row) = 0;
    end
end

disp('--- MISSION MAPPING RESULTS ---');
disp(sortrows(final_mapping, 'RepeatCorr', 'descend'));

% --- 4. Plotting ---
figure('Name', 'Final Verification', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);

for p = 1:numPatients
    subplot(4, 3, p);
    
    pID = final_mapping.PatientID(p);
    p_idx = find(strcmp(strrep(patientFolders, 'patient_', ''), char(pID)), 1);
    
    est_full = oae_results_cell{p_idx};
    t_ms = (0:length(est_full)-1)' / fs_values(p_idx) * 1000;
    
    crop_idx = (t_ms > 4 & t_ms < 16);
    est_crop = est_full(crop_idx);
    est_norm = est_crop / (max(abs(est_crop)) + eps);
    
    if final_mapping.Result(p) == "PASS" && final_mapping.Template(p) ~= "N/A"
        temp_raw = templates.(char(final_mapping.Template(p)))(:);
        
        temp_adj = [temp_raw; zeros(max(0, length(est_full) - length(temp_raw)), 1)];
        temp_adj = temp_adj(1:length(est_full));
        temp_crop = temp_adj(crop_idx);
        temp_norm = temp_crop / (max(abs(temp_crop)) + eps);
        
        [c_al, lags] = xcorr(est_norm, temp_norm, 'coeff');
        [~, b_lag] = max(c_al);
        temp_norm = circshift(temp_norm, lags(b_lag));
        
        plot(t_ms(crop_idx), est_norm, 'b', 'LineWidth', 1.3); hold on;
        plot(t_ms(crop_idx), temp_norm, 'r--', 'LineWidth', 1.3);
        title(sprintf('ID:%s | PASS | R:%.2f | S:%.1f', ...
            pID, final_mapping.RepeatCorr(p), final_mapping.Confidence(p)), 'FontSize', 8);
        legend('Estimated', 'Template', 'FontSize', 6, 'Location', 'best');
    else
        plot(t_ms(crop_idx), est_norm, 'k', 'LineWidth', 1.3);
        title(sprintf('ID:%s | REFER | R:%.2f', ...
            pID, final_mapping.RepeatCorr(p)), 'FontSize', 8);
    end
    
    grid on;
    xlim([4 16]);
    
    if p > 9
        xlabel('Time (ms)');
    end
    if mod(p,3) == 1
        ylabel('Norm Amp');
    end
end

sgtitle('Calculated OAE (Blue/Black) vs Reference Template (Red)', 'FontWeight', 'bold');