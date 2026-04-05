clear; clc;

% --- Configuration & Paths ---
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

fprintf('Processing %d patients... Please wait.\n', length(patientFolders));



for p = 1:length(patientFolders)
    pID = strrep(patientFolders{p}, 'patient_', '');
    pFolder = fullfile(basePath, patientFolders{p});
    
    % 1. Load Data
    [rec, fs] = audioread(fullfile(pFolder, ['Patient_' pID '_rec.wav']));
    info = jsondecode(fileread(fullfile(pFolder, ['Patient_' pID '_info.json'])));
    
    % 2. Extraction & Sub-Sample Alignment
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
                if isempty(p_ref); p_ref = seg(1:100); shifted = seg;
                else
                    [c_a, lags_a] = xcorr(seg(1:100), p_ref, 15, 'coeff');
                    [~, m_idx] = max(c_a);
                    shifted = circshift(seg, -lags_a(m_idx));
                end
                tmp = [tmp, shifted]; %#ok<AGROW>
            end
        end
        epochData{i} = tmp;
    end
    
    % 3. Nonlinear Summation
    minE = min(cellfun(@(x) size(x,2), epochData));
    acc = zeros(info.epochSize, 1);
    v_sets = 0;
    for k = 1:minE
        quad = epochData{1}(:,k) + epochData{2}(:,k) + epochData{3}(:,k) + epochData{4}(:,k);
        if all(isfinite(quad)); acc = acc + quad; v_sets = v_sets + 1; end
    end
    oae_raw = acc / max(1, v_sets);
    
    % 4. Wavelet-Boosted Post-Processing
    t = (0:info.epochSize-1)' / fs;
    
    % A. Initial Bandpass (Broad) to clear DC and high-hiss
    bpFilt = designfilt('bandpassiir','FilterOrder',6, ...
        'HalfPowerFrequency1',800,'HalfPowerFrequency2',4500, 'SampleRate',fs);
    pre_filtered = filtfilt(bpFilt, oae_raw);
    
    % B. Wavelet Denoising (The "Magic" Step)
    % 'db4' (Daubechies 4) is excellent for OAEs because its shape 
    % resembles a biological "wiggle".
    % 'UniversalThreshold' automatically calculates the noise floor.
    oae_denoised = wdenoise(pre_filtered, 5, ...
        'Wavelet', 'db4', ...
        'DenoisingMethod', 'UniversalThreshold', ...
        'ThresholdRule', 'Soft');
    
    % C. Apply Gaussian Window to focus on 4-15ms
    win_env = exp(-((t-0.009).^2)/(2*0.0045^2)); 
    oae_clean = oae_denoised .* win_env;

    % 5. Match against all templates
    m_idx = (t > 0.004 & t < 0.016);
    oae_crop = oae_clean(m_idx);
    for j = 1:length(temp_names)
        target = templates.(temp_names{j})(:);
        target_adj = [target; zeros(max(0, length(oae_clean)-length(target)), 1)];
        target_adj = target_adj(1:length(oae_clean));
        target_crop = target_adj(m_idx);
        [c, ~] = xcorr(oae_crop, target_crop, 'coeff');
        score = max(abs(c));
        if isnan(score); score = 0; end
        all_scores_matrix = [all_scores_matrix; p, j, score]; %#ok<AGROW>
    end
end



% ---  MAPPING LOGIC ---
all_scores_matrix = sortrows(all_scores_matrix, -3); 
final_mapping = table('Size', [length(patientFolders) 4], ...
    'VariableTypes', {'string', 'string', 'string', 'double'}, ...
    'VariableNames', {'PatientID', 'Result', 'Template', 'Confidence'});

for i = 1:length(patientFolders)
    final_mapping.PatientID(i) = string(strrep(patientFolders{i}, 'patient_', ''));
    final_mapping.Result(i) = "REFER";
    final_mapping.Template(i) = "N/A";
end

u_pat = []; u_temp = []; count = 0;
for i = 1:size(all_scores_matrix, 1)
    p_idx = all_scores_matrix(i, 1); t_idx = all_scores_matrix(i, 2); sc = all_scores_matrix(i, 3);
    if ~any(u_pat == p_idx) && ~any(u_temp == t_idx) && count < 8
        pID = string(strrep(patientFolders{p_idx}, 'patient_', ''));
        r_idx = find(final_mapping.PatientID == pID);
        final_mapping.Result(r_idx) = "PASS";
        final_mapping.Template(r_idx) = string(temp_names{t_idx});
        final_mapping.Confidence(r_idx) = sc * 100;
        u_pat = [u_pat; p_idx]; u_temp = [u_temp; t_idx]; count = count + 1;
    end
end


% confidence for REFERs
for i = 1:height(final_mapping)
    if final_mapping.Result(i) == "REFER"
        pID_str = final_mapping.PatientID(i);
        p_idx = find(strcmp(cellfun(@(x) strrep(x,'patient_',''), patientFolders, 'un', 0), pID_str));
        p_rows = all_scores_matrix(all_scores_matrix(:,1) == p_idx, :);
        final_mapping.Confidence(i) = p_rows(1,3) * 100;
    end
end

disp('--- FINAL MISSION MAPPING ---');
disp(sortrows(final_mapping, 'Confidence', 'descend'));

% --- 6. PLOT ALL ESTIMATED OAEs (Fixed Dimensions) ---
figure('Name', 'All Estimated OAEs', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for p = 1:length(patientFolders)
    subplot(4, 3, p);
    
    % Get the signal for this patient
    signal = oae_results_cell{p};
    
    % Re-calculate time vector for THIS specific signal length to avoid mismatch
    % (Using the fs from the last processed patient, or store fs in the loop)
    t_patient = (0:length(signal)-1)' / fs * 1000; 
    
    % Ensure t_patient and signal are the same length for plotting
    plot_len = min(length(t_patient), length(signal));
    
    plot(t_patient(1:plot_len), signal(1:plot_len), 'LineWidth', 1);
    
    % Get the result from our mapping table for the title
    pID_current = string(strrep(patientFolders{p}, 'patient_', ''));
    row_idx = find(final_mapping.PatientID == pID_current);
    
    res = char(final_mapping.Result(row_idx));
    conf = final_mapping.Confidence(row_idx);
    
    title(sprintf('ID: %s (%s)\nConf: %.1f%%', pID_current, res, conf), 'FontSize', 8);
    
    grid on; xlim([0 20]);
    if p > 9; xlabel('Time (ms)'); end
    if mod(p, 3) == 1; ylabel('Amp'); end
end
sgtitle('TEOAE Screening Results - All Patients', 'FontSize', 14, 'FontWeight', 'bold');