clear; clc; close all;

% --- 1. Configuration & Paths ---
basePath = '/Users/kourosh/Desktop/University/Self Study/Audio-Explorers20206/Diagnostics DSP/Patient Data';
templatePath = fullfile(basePath, 'lostOaes.json');

% Load Templates and Patient Folders 
templates = jsondecode(fileread(templatePath));
temp_names = fieldnames(templates);
dirInfo = dir(fullfile(basePath, 'patient_*'));
patientFolders = {dirInfo.name};
numPatients = length(patientFolders);

% Pre-allocate storage 
all_scores_matrix = []; % [PatientIdx, TemplateIdx, CombinedScore]
oae_results_cell = cell(numPatients, 1);
repeat_corr_values = zeros(numPatients, 1);
fs_values = zeros(numPatients, 1);

fprintf('Mission Start: Processing 12 newborns for DGS Titanius-12...\n'); 

% --- 2. Main Signal Recovery Loop ---
for p = 1:numPatients
    pID = strrep(patientFolders{p}, 'patient_', '');
    pFolder = fullfile(basePath, patientFolders{p});
    
    % A. Load Data & Metadata 
    [rec, fs] = audioread(fullfile(pFolder, ['Patient_' pID '_rec.wav']));
    fs_values(p) = fs;
    info = jsondecode(fileread(fullfile(pFolder, ['Patient_' pID '_info.json'])));
    
    % Get center frequencies for spectral validation 
    if isfield(info, 'centerFreq'), target_freqs = info.centerFreq(:);
    else, target_freqs = [1000; 2000; 3000; 4000]; end
    
    % B. Extraction & Sub-Sample Alignment (Handling Cosmic Jitter) 
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
                seg = rec(s:e) - mean(rec(s:e));
                % Correct for stochastic interference jitter 
                if isempty(p_ref), p_ref = seg(1:min(100, length(seg))); shifted = seg;
                else
                    [c_a, lags_a] = xcorr(seg(1:min(100, length(seg))), p_ref, 20, 'coeff');
                    [~, m_idx] = max(c_a);
                    shifted = circshift(seg, -lags_a(m_idx));
                end
                tmp = [tmp, shifted]; %#ok<AGROW>
            end
        end
        epochData{i} = tmp;
    end
    
    % C. Nonlinear Summation (+1+1+1-3 paradigm) 
    minE = min(cellfun(@(x) size(x,2), epochData));
    acc = zeros(info.epochSize, 1);
    acc1 = zeros(info.epochSize, 1); acc2 = zeros(info.epochSize, 1);
    v_sets = 0;
    
    for k = 1:minE
        quad = epochData{1}(:,k) + epochData{2}(:,k) + epochData{3}(:,k) + epochData{4}(:,k);
        if all(isfinite(quad))
            acc = acc + quad; v_sets = v_sets + 1;
            % Split-half for Reproducibility Check 
            if mod(k,2) == 1, acc1 = acc1 + quad; else, acc2 = acc2 + quad; end
        end
    end
    
    % D. Post-Processing (4-16ms Gating + Bandpass) 
    t = (0:info.epochSize-1)' / fs;
    win_env = exp(-((t - 0.008).^2) / (2 * 0.004^2)); 
    bpFilt = designfilt('bandpassiir', 'FilterOrder', 10, ...
        'HalfPowerFrequency1', 1200, 'HalfPowerFrequency2', 3800, 'SampleRate', fs);
    
    oae_clean = filtfilt(bpFilt, (acc/v_sets) .* win_env);
    oae_clean1 = filtfilt(bpFilt, (acc1/max(1,floor(v_sets/2))) .* win_env);
    oae_clean2 = filtfilt(bpFilt, (acc2/max(1,floor(v_sets/2))) .* win_env);
    
    % E. Noise Floor Compensation 
    noise_est = (mean(epochData{1},2) + mean(epochData{2},2)) - ...
                (mean(epochData{3},2) + mean(epochData{4},2));
    oae_clean = oae_clean - 0.1 * filtfilt(bpFilt, noise_est);
    oae_results_cell{p} = oae_clean;
    
    % F. Quality Metrics (Repeatability) 
    resp_idx = (t > 0.004 & t < 0.016);
    x1 = oae_clean1(resp_idx); x2 = oae_clean2(resp_idx);
    if norm(x1)>0 && norm(x2)>0, r_val = corr(x1, x2); else, r_val = 0; end
    repeat_corr_values(p) = max(0, r_val);
    
    % G. Dual-Domain Template Scoring (Phase-Sensitive) 
    oae_crop = oae_clean(resp_idx);
    Nfft = length(oae_crop);
    Ye_norm = abs(fft(oae_crop)/Nfft); Ye_norm = Ye_norm(1:floor(Nfft/2)+1)/max(Ye_norm);
    f_axis = fs * (0:floor(Nfft/2))/Nfft;
    
    for j = 1:length(temp_names)
        target = templates.(temp_names{j})(:);
        target_adj = [target; zeros(max(0, length(oae_clean)-length(target)), 1)];
        target_crop = target_adj(resp_idx);
        
        % 1. Time-Domain Morphology (No abs() to ensure phase match)
        [c, ~] = xcorr(oae_crop, target_crop, 'coeff');
        time_score = max(c); 
        
        % 2. Spectral Fit (Validating against CenterFreqs) 
        Yt_norm = abs(fft(target_crop)/Nfft); Yt_norm = Yt_norm(1:floor(Nfft/2)+1)/max(Yt_norm);
        spectral_fit = 0;
        for f_target = target_freqs'
            [~, bin] = min(abs(f_axis - f_target));
            spectral_fit = spectral_fit + (Ye_norm(bin) * Yt_norm(bin));
        end
        spectral_fit = spectral_fit / length(target_freqs);
        
        % 3. Combined Decision Logic
        all_scores_matrix = [all_scores_matrix; p, j, (0.7 * max(0,time_score) + 0.3 * spectral_fit)]; %#ok<AGROW>
    end
end

% --- 3. Competitive Global Assignment (Sudoku-Logic for Perfect Match)  ---
[~, rep_order] = sort(repeat_corr_values, 'descend');
pass_p_idx = rep_order(1:8); % Top 8 by reproducibility 

final_mapping = table('Size', [numPatients 5], ...
    'VariableTypes', {'string', 'string', 'string', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'Result', 'Template', 'Confidence', 'RepeatCorr'});

for i = 1:numPatients
    final_mapping.PatientID(i) = string(strrep(patientFolders{i}, 'patient_', ''));
    final_mapping.Result(i) = "REFER"; final_mapping.Template(i) = "N/A";
    final_mapping.RepeatCorr(i) = repeat_corr_values(i);
end

% Build Cost Matrix [8 PASS Patients x 8 Templates]
cost_matrix = zeros(8, 8);
for i = 1:8
    p_orig = pass_p_idx(i);
    p_scores = all_scores_matrix(all_scores_matrix(:,1) == p_orig, 3);
    cost_matrix(i, :) = p_scores';
end

% Hungarian-Style Global Best Match 
temp_cost = cost_matrix;
for k = 1:8
    [best_val, lin_idx] = max(temp_cost(:));
    [r, c] = ind2sub(size(temp_cost), lin_idx);
    
    p_row = find(final_mapping.PatientID == string(strrep(patientFolders{pass_p_idx(r)}, 'patient_', '')));
    final_mapping.Result(p_row) = "PASS"; 
    final_mapping.Template(p_row) = string(temp_names{c}); 
    final_mapping.Confidence(p_row) = best_val * 100;
    
    temp_cost(r, :) = -inf; % Competitive lockout
    temp_cost(:, c) = -inf;
end

% 4. Final Output & Morphology Check
disp('--- FINAL MISSION RESULTS: DGS TITANIUS-12 ---'); 
disp(sortrows(final_mapping, 'RepeatCorr', 'descend'));

figure('Name', 'Final Validation: Blue (OAE) vs Red (Template)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
for p = 1:numPatients
    subplot(4, 3, p);
    pID = final_mapping.PatientID(p);
    p_orig = find(contains(patientFolders, pID));
    sig = oae_results_cell{p_orig};
    t_plot = (0:length(sig)-1)' / fs_values(p_orig) * 1000;
    win = (t_plot > 4 & t_plot < 16);
    
    s_norm = sig(win) / (max(abs(sig(win))) + eps);
    if final_mapping.Result(p) == "PASS"
        tmp = templates.(char(final_mapping.Template(p)))(:);
        tmp_crop = tmp(1:min(length(tmp), sum(win)));
        t_norm = tmp_crop / (max(abs(tmp_crop)) + eps);
        [c_al, lags] = xcorr(s_norm, t_norm, 'coeff'); [~, bl] = max(c_al);
        plot(t_plot(win), s_norm, 'b'); hold on;
        plot(t_plot(win), circshift(t_norm, lags(bl)), 'r--');
        title(sprintf('ID:%s | R:%.2f', pID, final_mapping.RepeatCorr(p)));
    else
        plot(t_plot(win), s_norm, 'k'); title(sprintf('ID:%s (REFER)', pID));
    end
    grid on; xlim([4 16]);
end
sgtitle('Recovered OAEs vs Shuffled Templates', 'FontWeight', 'bold');


% --- 5. EXPORT FOR GUI ---
exportDir = fullfile(basePath, 'gui_export');
if ~exist(exportDir, 'dir')
    mkdir(exportDir);
end

% Save final mapping
writetable(final_mapping, fullfile(exportDir, 'final_mapping.csv'));

% Save template scores
scores_table = array2table(all_scores_matrix, ...
    'VariableNames', {'PatientIdx','TemplateIdx','Score'});

scores_table.PatientID = strings(height(scores_table),1);
scores_table.Template = strings(height(scores_table),1);

for i = 1:height(scores_table)
    scores_table.PatientID(i) = string(strrep(patientFolders{scores_table.PatientIdx(i)}, 'patient_', ''));
    scores_table.Template(i) = string(temp_names{scores_table.TemplateIdx(i)});
end

scores_table = scores_table(:, {'PatientID','Template','Score'});
writetable(scores_table, fullfile(exportDir, 'template_scores.csv'));

% Build patient_results struct for GUI
patient_results = struct([]);

for p = 1:numPatients
    pID = string(strrep(patientFolders{p}, 'patient_', ''));
    row_idx = find(final_mapping.PatientID == pID, 1);

    sig = oae_results_cell{p};
    t_ms = (0:length(sig)-1)' / fs_values(p) * 1000;
    crop_idx = (t_ms > 4 & t_ms < 16);

    est_crop = sig(crop_idx);
    est_norm = est_crop / (max(abs(est_crop)) + eps);

    assigned_template = char(final_mapping.Template(row_idx));
    matched_norm = zeros(size(est_norm));
    template_fft_freq = 0;
    template_fft_mag = 0;

    if final_mapping.Result(row_idx) == "PASS" && final_mapping.Template(row_idx) ~= "N/A"
        tmp = templates.(assigned_template)(:);
        tmp_adj = [tmp; zeros(max(0, length(sig)-length(tmp)), 1)];
        tmp_adj = tmp_adj(1:length(sig));
        tmp_crop = tmp_adj(crop_idx);
        tmp_norm = tmp_crop / (max(abs(tmp_crop)) + eps);

        [c_al, lags] = xcorr(est_norm, tmp_norm, 'coeff');
        [~, bl] = max(c_al);
        matched_norm = circshift(tmp_norm, lags(bl));

        Nfft_t = length(tmp_crop);
        Yt = abs(fft(tmp_crop) / Nfft_t);
        template_fft_mag = Yt(1:floor(Nfft_t/2)+1) / (max(Yt) + eps);
        template_fft_freq = fs_values(p) * (0:floor(Nfft_t/2)) / Nfft_t;
    end

    Nfft = length(est_crop);
    Ye = abs(fft(est_crop) / Nfft);
    fft_mag = Ye(1:floor(Nfft/2)+1) / (max(Ye) + eps);
    fft_freq = fs_values(p) * (0:floor(Nfft/2)) / Nfft;

    patient_results(p).PatientID = char(pID);
    patient_results(p).fs = fs_values(p);
    patient_results(p).t_ms = t_ms;
    patient_results(p).t_crop_ms = t_ms(crop_idx);
    patient_results(p).oae_clean = sig;
    patient_results(p).oae_crop = est_crop;
    patient_results(p).est_norm = est_norm;
    patient_results(p).matched_norm = matched_norm;
    patient_results(p).fft_freq = fft_freq;
    patient_results(p).fft_mag = fft_mag;
    patient_results(p).template_fft_freq = template_fft_freq;
    patient_results(p).template_fft_mag = template_fft_mag;
    patient_results(p).result = char(final_mapping.Result(row_idx));
    patient_results(p).assigned_template = assigned_template;
    patient_results(p).assigned_score = final_mapping.Confidence(row_idx);
    patient_results(p).repeat_corr = final_mapping.RepeatCorr(row_idx);
end

save(fullfile(exportDir, 'patient_results.mat'), 'patient_results', 'temp_names');

% Save summary JSON
summary.final_mapping_file = 'final_mapping.csv';
summary.template_scores_file = 'template_scores.csv';
summary.patient_results_file = 'patient_results.mat';
summary.num_patients = numPatients;
summary.template_names = temp_names;

jsonText = jsonencode(summary, PrettyPrint=true);
fid = fopen(fullfile(exportDir, 'summary.json'), 'w');
fprintf(fid, '%s', jsonText);
fclose(fid);

disp('--- GUI EXPORT COMPLETE ---');
disp(exportDir);