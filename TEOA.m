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
    
% --- 4. Boosted Post-Processing (Time-Frequency Gating) ---
    t = (0:info.epochSize-1)' / fs;
    
    % Create a Time-Frequency Gate (Narrower at the start, wider at the end)
    % This kills early low-freq noise and late high-freq noise
    oae_boosted = oae_raw;
    win_env = exp(-((t-0.008).^2)/(2*0.004^2)); % Gaussian centered at 8ms
    oae_boosted = oae_boosted .* win_env;
    
    % Zero-phase filtering with a tighter transition
    bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
        'HalfPowerFrequency1',1200,'HalfPowerFrequency2',3800, 'SampleRate',fs);
    oae_clean = filtfilt(bpFilt, oae_boosted);
    
    % Final Noise Floor Leveling
    noise_est = (mean(epochData{1},2) + mean(epochData{2},2)) - ...
                (mean(epochData{3},2) + mean(epochData{4},2));
    oae_clean = oae_clean - (0.1 * filtfilt(bpFilt, noise_est)); % Subtract 10% of estimated noise
    
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










t = (0:epochSize-1)/fs*1000; % ms

figure;
plot(t, oae_est);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['Estimated OAE for patient ' patient]);
grid on;

oae_est = oae_est(:);

N = length(oae_est);
f = (0:N-1)*(fs/N);

Y = fft(oae_est);
magY = abs(Y);

figure;
plot(f(1:floor(N/2)), magY(1:floor(N/2)), 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('|FFT|');
title('Fourier Transform of Estimated OAE');



templates = jsondecode(fileread('lostOaes.json'));
names = fieldnames(templates);
% pick one template (e.g., first one)
temp = templates.(names{6});
temp = temp(:);
% FFT
N = length(temp);
Y = fft(temp);
P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:floor(N/2))/N;
figure;
plot(f, P1, 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(['FFT of template: ' names{1}]);






% figure;
% hold on;
% 
% for k = 1:min(10, size(alignedEpochs{1},2))
%     x = alignedEpochs{1}(:,k);
%     x = x - mean(x);
% 
%     N = length(x);
%     Y = fft(x);
%     P2 = abs(Y/N);
%     P1 = P2(1:floor(N/2)+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = fs*(0:floor(N/2))/N;
% 
%     plot(f, P1);
% end
% 
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('FFT of individual aligned epochs (Type A)');
% xlim([0 5000]);

% figure;
% hold on;
% 
% 
% for k = 1:min(10, size(alignedEpochs{2},2))
%     x = alignedEpochs{2}(:,k);
%     x = x - mean(x);
% 
%     N = length(x);
%     Y = fft(x);
%     P2 = abs(Y/N);
%     P1 = P2(1:floor(N/2)+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = fs*(0:floor(N/2))/N;
% 
%     plot(f, P1);
% end
% 
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('FFT of individual aligned epochs (Type A)');
% xlim([0 5000]);







% 
% templates = jsondecode(fileread('lostOaes.json'));
% names = fieldnames(templates);
% 
% W = [1 1 1 -3];
% 
% alignStartList = [20 30 40 50 60];
% alignEndList   = [120 160 200 240 280];
% maxLagList     = [10 20 30 40 50 60 70 80];
% 
% bestScore = -inf;
% bestParams = struct();
% bestOAE = [];
% bestTemplateName = '';
% bestTemplateScore = -inf;
% bestRepeatability = -inf;
% 
% for alignStart = alignStartList
%     for alignEnd = alignEndList
%         if alignEnd <= alignStart || alignEnd > epochSize
%             continue;
%         end
% 
%         for maxLag = maxLagList
% 
%             % -------- Align epochs and compute avgEpoch --------
%             avgEpoch = cell(1,4);
%             alignedEpochs = cell(1,4);
% 
%             ok = true;
% 
%             for i = 1:4
%                 X = epochs{i};
% 
%                 if isempty(X) || size(X,2) < 4
%                     ok = false;
%                     break;
%                 end
% 
%                 ref = median(X, 2);
% 
%                 for iter = 1:3
%                     X_aligned = zeros(size(X));
%                     refWin = ref(alignStart:alignEnd);
% 
%                     for k = 1:size(X,2)
%                         x = X(:,k);
%                         xWin = x(alignStart:alignEnd);
% 
%                         [xc, lags] = xcorr(xWin, refWin, maxLag, 'coeff');
%                         [~, idxMax] = max(abs(xc));
%                         lag = lags(idxMax);
% 
%                         if lag > 0
%                             X_aligned(:,k) = [x(lag+1:end); zeros(lag,1)];
%                         elseif lag < 0
%                             s = -lag;
%                             X_aligned(:,k) = [zeros(s,1); x(1:end-s)];
%                         else
%                             X_aligned(:,k) = x;
%                         end
%                     end
% 
%                     X = X_aligned;
%                     ref = median(X, 2);
%                 end
% 
%                 alignedEpochs{i} = X_aligned;
%                 avgEpoch{i} = median(X_aligned, 2);
%             end
% 
%             if ~ok
%                 continue;
%             end
% 
%             % -------- Fixed nonlinear combination --------
%             oae_est = W(1)*avgEpoch{1} + W(2)*avgEpoch{2} + ...
%                       W(3)*avgEpoch{3} + W(4)*avgEpoch{4};
% 
%             oae_est = oae_est(:);
%             oae_est = oae_est - mean(oae_est);
% 
%             if norm(oae_est) == 0
%                 continue;
%             end
% 
%             x = oae_est / norm(oae_est);
% 
%             % -------- Template score --------
%             localBestTemplateScore = -inf;
%             localBestTemplateName = '';
% 
%             for k = 1:length(names)
%                 temp = templates.(names{k});
%                 temp = temp(:);
% 
%                 temp_rs = resample(temp, length(x), length(temp));
% 
%                 if norm(temp_rs) == 0
%                     continue;
%                 end
% 
%                 temp_rs = temp_rs / norm(temp_rs);
% 
%                 [xc, ~] = xcorr(x, temp_rs, 'coeff');
%                 s = max(abs(xc));
% 
%                 if s > localBestTemplateScore
%                     localBestTemplateScore = s;
%                     localBestTemplateName = names{k};
%                 end
%             end
% 
%             % -------- Odd/even repeatability --------
%             oddAvg = cell(1,4);
%             evenAvg = cell(1,4);
% 
%             validRep = true;
%             for i = 1:4
%                 X = alignedEpochs{i};
% 
%                 if size(X,2) < 4
%                     validRep = false;
%                     break;
%                 end
% 
%                 oddAvg{i}  = median(X(:,1:2:end), 2);
%                 evenAvg{i} = median(X(:,2:2:end), 2);
%             end
% 
%             if ~validRep
%                 continue;
%             end
% 
%             oae_odd = W(1)*oddAvg{1} + W(2)*oddAvg{2} + ...
%                       W(3)*oddAvg{3} + W(4)*oddAvg{4};
% 
%             oae_even = W(1)*evenAvg{1} + W(2)*evenAvg{2} + ...
%                        W(3)*evenAvg{3} + W(4)*evenAvg{4};
% 
%             oae_odd = oae_odd(:) - mean(oae_odd);
%             oae_even = oae_even(:) - mean(oae_even);
% 
%             if norm(oae_odd) == 0 || norm(oae_even) == 0
%                 continue;
%             end
% 
%             xo = oae_odd / norm(oae_odd);
%             xe = oae_even / norm(oae_even);
% 
%             [xcRep, ~] = xcorr(xo, xe, 'coeff');
%             repeatability = max(abs(xcRep));
% 
%             % -------- Final score --------
%             score = 0.7 * localBestTemplateScore + 0.3 * repeatability;
% 
%             if score > bestScore
%                 bestScore = score;
%                 bestParams.alignStart = alignStart;
%                 bestParams.alignEnd = alignEnd;
%                 bestParams.maxLag = maxLag;
% 
%                 bestOAE = x;
%                 bestTemplateName = localBestTemplateName;
%                 bestTemplateScore = localBestTemplateScore;
%                 bestRepeatability = repeatability;
%             end
%         end
%     end
% end
% 
% disp('Best parameters found:')
% disp(bestParams)
% fprintf('Best overall score      = %.4f\n', bestScore);
% fprintf('Best template score     = %.4f\n', bestTemplateScore);
% fprintf('Best repeatability      = %.4f\n', bestRepeatability);
% fprintf('Best template candidate = %s\n', bestTemplateName);




for i = 1:4
    idx = find(trigger{i} ~= 0);
    fprintf('\nType %s\n', types{i});
    fprintf('Unique nonzero trigger values:\n');
    disp(unique(trigger{i}(idx))');
end