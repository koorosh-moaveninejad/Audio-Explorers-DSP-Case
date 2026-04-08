clear; clc;

% Configuration & Paths
basePath = '/Users/kourosh/Desktop/University/Self Study/Audio-Explorers20206/Diagnostics DSP/Patient Data';
templatePath = fullfile(basePath, 'lostOaes.json');

% Load Templates and Folders
templates = jsondecode(fileread(templatePath));
temp_names = fieldnames(templates);
dirInfo = dir(fullfile(basePath, 'patient_*'));
patientFolders = {dirInfo.name};

% Pre-allocate storage
all_scores_matrix = []; % [PatientIdx, TemplateIdx, Score]
oae_results_cell = cell(length(patientFolders), 1); % Store waveforms for plotting
fft_results_cell = cell(length(patientFolders), 1);
fft_template_cell = cell(length(patientFolders), 1);
freq_axis_cell = cell(length(patientFolders), 1);
snr_values = zeros(length(patientFolders),1);     % SNR per patient
repeat_corr_values = zeros(length(patientFolders),1); % for split-half reproducibility
fs_values = zeros(length(patientFolders),1);      
fprintf('Processing %d patients... Please wait.\n', length(patientFolders));


for p = 1:length(patientFolders)
    pID = strrep(patientFolders{p}, 'patient_', '');
    pFolder = fullfile(basePath, patientFolders{p});
    
    % 1. Load Data
    [rec, fs] = audioread(fullfile(pFolder, ['Patient_' pID '_rec.wav']));
    fs_values(p) = fs;
    info = jsondecode(fileread(fullfile(pFolder, ['Patient_' pID '_info.json'])));
    
    % 2. Extraction & Alignment
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
                
                if isempty(p_ref)
                    p_ref = seg(1:min(100,length(seg)));
                    shifted = seg;
                else
                    ref_part = p_ref(1:min(length(p_ref),100));
                    seg_part = seg(1:min(length(seg),length(ref_part)));
                    
                    [c_a, lags_a] = xcorr(seg_part, ref_part, 15, 'coeff');
                    [~, m_idx_align] = max(c_a);
                    shifted = circshift(seg, -lags_a(m_idx_align));
                end
                
                tmp = [tmp, shifted]; %#ok<AGROW>
            end
        end
        
        epochData{i} = tmp;
    end
    
    % 3. Nonlinear Summation + split-half reproducibility
    minE = min(cellfun(@(x) size(x,2), epochData));
    
    acc  = zeros(info.epochSize, 1); v_sets  = 0;
    acc1 = zeros(info.epochSize, 1); v_sets1 = 0;
    acc2 = zeros(info.epochSize, 1); v_sets2 = 0;
    
    for k = 1:minE
        quad = epochData{1}(:,k) + epochData{2}(:,k) + epochData{3}(:,k) + epochData{4}(:,k);
        
        if all(isfinite(quad))
            acc = acc + quad;
            v_sets = v_sets + 1;
            
            if mod(k,2) == 1
                acc1 = acc1 + quad;
                v_sets1 = v_sets1 + 1;
            else
                acc2 = acc2 + quad;
                v_sets2 = v_sets2 + 1;
            end
        end
    end
    
    oae_raw  = acc  / max(1, v_sets);
    oae_raw1 = acc1 / max(1, v_sets1);
    oae_raw2 = acc2 / max(1, v_sets2);
    
    % 4. Post-Processing
    t = (0:info.epochSize-1)' / fs;
    
    % Time gate
    win_env = exp(-((t-0.008).^2) / (2*0.004^2));  % centered at 8 ms
    
    % Bandpass
    bpFilt = designfilt('bandpassiir', 'FilterOrder', 10, ...
        'HalfPowerFrequency1', 1200, ...
        'HalfPowerFrequency2', 3800, ...
        'SampleRate', fs);
    
    % Full estimate
    oae_boosted = oae_raw .* win_env;
    oae_clean = filtfilt(bpFilt, oae_boosted);
    
    % Split-half estimate 1
    oae_boosted1 = oae_raw1 .* win_env;
    oae_clean1 = filtfilt(bpFilt, oae_boosted1);
    
    % Split-half estimate 2
    oae_boosted2 = oae_raw2 .* win_env;
    oae_clean2 = filtfilt(bpFilt, oae_boosted2);
    
    % Noise estimate
    noise_est = (mean(epochData{1},2) + mean(epochData{2},2)) - ...
                (mean(epochData{3},2) + mean(epochData{4},2));
    noise_clean = filtfilt(bpFilt, noise_est);
    
    % Final leveling
    oae_clean  = oae_clean  - (0.1 * noise_clean);
    oae_clean1 = oae_clean1 - (0.1 * noise_clean);
    oae_clean2 = oae_clean2 - (0.1 * noise_clean);
    
    oae_results_cell{p} = oae_clean;
    % --- Fourier analysis of estimated OAE in response window ---
    resp_idx = (t > 0.004 & t < 0.016);
    oae_fft_seg = oae_clean(resp_idx);
    
    Nfft = length(oae_fft_seg);
    if Nfft > 1
        Y = fft(oae_fft_seg);
        P2 = abs(Y / Nfft);
        P1 = P2(1:floor(Nfft/2)+1);
        if length(P1) > 2
            P1(2:end-1) = 2 * P1(2:end-1);
        end
        f_axis = fs * (0:floor(Nfft/2)) / Nfft;
    else
        P1 = 0;
        f_axis = 0;
    end
    
    fft_results_cell{p} = P1;
    freq_axis_cell{p} = f_axis;
    
    % Zero-phase filtering with a tighter transition
    bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
        'HalfPowerFrequency1',1200,'HalfPowerFrequency2',3800, 'SampleRate',fs);
    oae_clean = filtfilt(bpFilt, oae_boosted);
    
    % Final Noise Floor Leveling
    noise_est = (mean(epochData{1},2) + mean(epochData{2},2)) - ...
                (mean(epochData{3},2) + mean(epochData{4},2));
    oae_clean = oae_clean - (0.1 * filtfilt(bpFilt, noise_est)); % Subtract 10% of estimated noise
    oae_results_cell{p} = oae_clean;

    % 5. Quality Metrics: SNR + reproducibility
    resp_idx  = (t > 0.004 & t < 0.016);   % response window
    noise_idx = (t > 0.017 & t < 0.020);   % late noise window
    
    resp_rms = rms(oae_clean(resp_idx));
    noise_rms = rms(oae_clean(noise_idx));
    snr_db = 20 * log10(resp_rms / max(noise_rms, eps));
    snr_values(p) = snr_db;
    
    x1 = oae_clean1(resp_idx);
    x2 = oae_clean2(resp_idx);
    
    if norm(x1) > 0 && norm(x2) > 0
        repeat_corr = corr(x1, x2);
    else
        repeat_corr = 0;
    end
    
    if isnan(repeat_corr)
        repeat_corr = 0;
    end
    repeat_corr_values(p) = repeat_corr;


    % 6. Match against all templates
    oae_crop = oae_clean(resp_idx);
    
    for j = 1:length(temp_names)
        target = templates.(temp_names{j})(:);
        
        target_adj = [target; zeros(max(0, length(oae_clean)-length(target)), 1)];
        target_adj = target_adj(1:length(oae_clean));
        target_crop = target_adj(resp_idx);
        % FFT of matched template crop
        Nfft_t = length(target_crop);
        if Nfft_t > 1
            Yt = fft(target_crop);
            P2t = abs(Yt / Nfft_t);
            P1t = P2t(1:floor(Nfft_t/2)+1);
            if length(P1t) > 2
                P1t(2:end-1) = 2 * P1t(2:end-1);
            end
            f_axis_t = fs_values(p) * (0:floor(Nfft_t/2)) / Nfft_t;
        else
            P1t = 0;
            f_axis_t = 0;
        end
        
        fft_template_cell{p} = P1t;
        
        [c, ~] = xcorr(oae_crop, target_crop, 'coeff');
        score = max(c);
        
        if isnan(score) || score < 0
            score = 0;
        end
        
        all_scores_matrix = [all_scores_matrix; p, j, score]; %#ok<AGROW>
    end
end



% 7. MAPPING LOGIC BASED ON REPEAT CORRELATION
all_scores_matrix = sortrows(all_scores_matrix, -3);

final_mapping = table('Size', [length(patientFolders) 5], ...
    'VariableTypes', {'string', 'string', 'string', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'Result', 'Template', 'Confidence', 'RepeatCorr'});

for i = 1:length(patientFolders)
    final_mapping.PatientID(i) = string(strrep(patientFolders{i}, 'patient_', ''));
    final_mapping.Result(i) = "REFER";
    final_mapping.Template(i) = "N/A";
    final_mapping.Confidence(i) = 0;
    final_mapping.RepeatCorr(i) = repeat_corr_values(i);
end

% Step 1: choose PASS patients by RepeatCorr
[~, rep_order] = sort(repeat_corr_values, 'descend');
pass_count = min(8, length(patientFolders));   % choose top 8 by RepeatCorr
pass_patients = rep_order(1:pass_count);

% Step 2: assign templates only among PASS patients
used_templates = [];

for n = 1:length(pass_patients)
    p_idx = pass_patients(n);

    % get all template matches for this patient
    p_rows = all_scores_matrix(all_scores_matrix(:,1) == p_idx, :);

    assigned = false;
    for r = 1:size(p_rows,1)
        t_idx = p_rows(r,2);
        sc    = p_rows(r,3);

        if ~any(used_templates == t_idx)
            pID = string(strrep(patientFolders{p_idx}, 'patient_', ''));
            row_idx = find(final_mapping.PatientID == pID);

            final_mapping.Result(row_idx) = "PASS";
            final_mapping.Template(row_idx) = string(temp_names{t_idx});
            final_mapping.Confidence(row_idx) = sc * 100;

            used_templates = [used_templates; t_idx]; %#ok<AGROW>
            assigned = true;
            break;
        end
    end

    % if all templates already used, still keep as PASS but no template
    if ~assigned
        pID = string(strrep(patientFolders{p_idx}, 'patient_', ''));
        row_idx = find(final_mapping.PatientID == pID);

        final_mapping.Result(row_idx) = "PASS";
        final_mapping.Template(row_idx) = "N/A";

        if ~isempty(p_rows)
            final_mapping.Confidence(row_idx) = p_rows(1,3) * 100;
        end
    end
end

% Step 3: REFER patients keep their best template score only as confidence
for i = 1:height(final_mapping)
    if final_mapping.Result(i) == "REFER"
        pID_str = final_mapping.PatientID(i);
        p_idx = find(strcmp(cellfun(@(x) strrep(x,'patient_',''), patientFolders, 'UniformOutput', false), pID_str));
        p_rows = all_scores_matrix(all_scores_matrix(:,1) == p_idx, :);

        if ~isempty(p_rows)
            final_mapping.Confidence(i) = p_rows(1,3) * 100;
        end
    end
end

% OAE existence rule 

[~, rep_order] = sort(repeat_corr_values, 'descend');
oae_exists_rule = false(length(repeat_corr_values),1);
oae_exists_rule(rep_order(1:min(8,end))) = true;

disp('--- TEMPLATE-INDEPENDENT OAE EXISTENCE CHECK ---');
existence_table = table( ...
    string(cellfun(@(x) strrep(x,'patient_',''), patientFolders, 'UniformOutput', false))', ...
    repeat_corr_values, oae_exists_rule, ...
    'VariableNames', {'PatientID','RepeatCorr','OAE_Exists'});

disp(sortrows(existence_table, 'RepeatCorr', 'descend'));

% 8. PLOT ALL ESTIMATED OAEs 
figure('Name', 'All Estimated OAEs', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for p = 1:length(patientFolders)
    subplot(4, 3, p);
    
    signal = oae_results_cell{p};
    t_patient = (0:length(signal)-1)' / fs_values(p) * 1000;
    
    plot_len = min(length(t_patient), length(signal));
    plot(t_patient(1:plot_len), signal(1:plot_len), 'LineWidth', 1);
    
    pID_current = string(strrep(patientFolders{p}, 'patient_', ''));
    row_idx = find(final_mapping.PatientID == pID_current);
    
    res = char(final_mapping.Result(row_idx));
    conf = final_mapping.Confidence(row_idx);
    % snr_here = final_mapping.SNR_dB(row_idx);
    rep_here = final_mapping.RepeatCorr(row_idx);
    
    title(sprintf('ID: %s (%s)\nConf: %.1f%% | R: %.2f', ...
        pID_current, res, conf, rep_here), 'FontSize', 8);
    
    grid on;
    xlim([0 20]);
    
    if p > 9
        xlabel('Time (ms)');
    end
    if mod(p,3) == 1
        ylabel('Amp');
    end
end

sgtitle('TEOAE Screening Results - All Patients', 'FontSize', 14, 'FontWeight', 'bold');

% 9. PLOT ESTIMATED OAE VS MATCHED TEMPLATE 
figure('Name', 'Estimated vs Matched Template OAEs', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);

for p = 1:length(patientFolders)
    subplot(4, 3, p);
    
    pID_current = string(strrep(patientFolders{p}, 'patient_', ''));
    row_idx = find(final_mapping.PatientID == pID_current);
    
    est_full = oae_results_cell{p};
    t_est = (0:length(est_full)-1)' / fs_values(p);
    
    m_idx_est = (t_est > 0.004 & t_est < 0.016);
    est_crop = est_full(m_idx_est);
    
    if norm(est_crop) > 0
        est_norm = est_crop / norm(est_crop);
    else
        est_norm = est_crop;
    end
    
    matched_norm = zeros(size(est_norm));
    template_name = char(final_mapping.Template(row_idx));
    
    if final_mapping.Template(row_idx) ~= "N/A"
        target = templates.(template_name)(:);
        
        target_adj = [target; zeros(max(0, length(est_full)-length(target)), 1)];
        target_adj = target_adj(1:length(est_full));
        target_crop = target_adj(m_idx_est);
        
        if norm(target_crop) > 0
            target_norm = target_crop / norm(target_crop);
        else
            target_norm = target_crop;
        end
        
        [c_align, lags_align] = xcorr(est_norm, target_norm, 'coeff');
        [~, best_idx] = max(c_align);
        best_lag = lags_align(best_idx);
        
        matched_norm = circshift(target_norm, best_lag);
    end
    
    t_crop_ms = t_est(m_idx_est) * 1000;
    
    plot(t_crop_ms, est_norm, 'LineWidth', 1.5); hold on;
    plot(t_crop_ms, matched_norm, '--', 'LineWidth', 1.5);
    hold off;
    
    grid on;
    xlim([4 16]);
    
    res = char(final_mapping.Result(row_idx));
    conf = final_mapping.Confidence(row_idx);
    % snr_here = final_mapping.SNR_dB(row_idx);
    rep_here = final_mapping.RepeatCorr(row_idx);
    
    title(sprintf('ID: %s (%s)\n%s | Conf: %.1f%% | R: %.2f', ...
        pID_current, res, template_name, conf, rep_here), 'FontSize', 8);
    
    if p > 9
        xlabel('Time (ms)');
    end
    if mod(p,3) == 1
        ylabel('Normalized Amp');
    end
    
    legend('Estimated', 'Matched Template', 'FontSize', 7, 'Location', 'best');
end

sgtitle('Estimated OAEs vs Matched Templates', 'FontSize', 14, 'FontWeight', 'bold');

 % --- 10. PLOT SNR VS REPRODUCIBILITY ---
% figure('Name', 'SNR vs Reproducibility', 'Color', 'w');
% 
% scatter(snr_values, repeat_corr_values, 80, 'filled');
% grid on;
% xlabel('SNR (dB)');
% ylabel('Repeatability Correlation');

% for p = 1:length(patientFolders)
%     pID_current = strrep(patientFolders{p}, 'patient_', '');
%     text(repeat_corr_values(p), pID_current, 'FontSize', 8);
% end
% 
% title('Patient Quality Metrics: SNR vs Split-Half Reproducibility');
% xline(3, '--');
% yline(0.35, '--');


% --- 11. PLOT FFT OF ESTIMATED OAE VS MATCHED TEMPLATE ---
figure('Name', 'FFT: Estimated vs Matched Template', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);

for p = 1:length(patientFolders)
    subplot(4, 3, p);

    pID_current = string(strrep(patientFolders{p}, 'patient_', ''));
    row_idx = find(final_mapping.PatientID == pID_current);

    % Estimated OAE FFT
    est_full = oae_results_cell{p};
    t_est = (0:length(est_full)-1)' / fs_values(p);
    resp_idx = (t_est > 0.004 & t_est < 0.016);
    est_crop = est_full(resp_idx);

    Nfft = length(est_crop);
    if Nfft > 1
        Ye = fft(est_crop);
        P2e = abs(Ye / Nfft);
        P1e = P2e(1:floor(Nfft/2)+1);
        if length(P1e) > 2
            P1e(2:end-1) = 2 * P1e(2:end-1);
        end
        f_est = fs_values(p) * (0:floor(Nfft/2)) / Nfft;
    else
        P1e = 0;
        f_est = 0;
    end

    % Default template FFT
    P1t = zeros(size(P1e));
    f_temp = f_est;

    template_name = char(final_mapping.Template(row_idx));
    if final_mapping.Template(row_idx) ~= "N/A"
        target = templates.(template_name)(:);
        target_adj = [target; zeros(max(0, length(est_full)-length(target)), 1)];
        target_adj = target_adj(1:length(est_full));
        target_crop = target_adj(resp_idx);

        Nfft_t = length(target_crop);
        if Nfft_t > 1
            Yt = fft(target_crop);
            P2t = abs(Yt / Nfft_t);
            P1t = P2t(1:floor(Nfft_t/2)+1);
            if length(P1t) > 2
                P1t(2:end-1) = 2 * P1t(2:end-1);
            end
            f_temp = fs_values(p) * (0:floor(Nfft_t/2)) / Nfft_t;
        end
    end

    plot(f_est, P1e, 'LineWidth', 1.5); hold on;
    plot(f_temp, P1t, '--', 'LineWidth', 1.5); hold off;

    grid on;
    xlim([500 5000]);

    res = char(final_mapping.Result(row_idx));
    conf = final_mapping.Confidence(row_idx);
    rep_here = final_mapping.RepeatCorr(row_idx);

    title(sprintf('ID: %s (%s)\n%s | Conf: %.1f%% | R: %.2f', ...
        pID_current, res, template_name, conf, rep_here), 'FontSize', 8);

    if p > 9
        xlabel('Frequency (Hz)');
    end
    if mod(p,3) == 1
        ylabel('Magnitude');
    end

    legend('Estimated FFT', 'Template FFT', 'FontSize', 7, 'Location', 'best');
end

sgtitle('FFT Comparison: Estimated OAE vs Matched Template', ...
    'FontSize', 14, 'FontWeight', 'bold');