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



for i = 1:4
    trigger{i} = audioread(fullfile(patientFolder, ['Patient_' patient '_trigger' types{i} '.wav']));
    active{i}  = audioread(fullfile(patientFolder, ['Patient_' patient '_active'  types{i} '.wav']));
end

info = jsondecode(fileread(fullfile(patientFolder, ['Patient_' patient '_info.json'])));
epochSize = info.epochSize;
rec = rec / info.recNorm;







activeStarts = cell(1,4);
activeEnds   = cell(1,4);
epochSigns   = cell(1,4);
epochs       = cell(1,4);

for i = 1:4
    actMask = active{i} > 0.5;

    % active regions
    aStarts = find(diff([0; actMask]) == 1);
    aEnds   = find(diff([actMask; 0]) == -1);

    activeStarts{i} = aStarts;
    activeEnds{i}   = aEnds;

    % all nonzero trigger events with sign
    trigIdxAll = find(trigger{i} ~= 0);
    trigValAll = trigger{i}(trigIdxAll);

    nBlocks = length(aStarts);
    tmpEpochs = zeros(epochSize, nBlocks);
    tmpSigns  = zeros(nBlocks,1);
    validCount = 0;

    for k = 1:nBlocks
        s = aStarts(k);
        e = s + epochSize - 1;

        if e > length(rec)
            continue;
        end

        % trigger events inside this active block
        inside = trigIdxAll >= aStarts(k) & trigIdxAll <= aEnds(k);

        if ~any(inside)
            continue;
        end

        % usually there should be one trigger in the block;
        % if more than one, take the first one
        trigPos = trigIdxAll(find(inside, 1, 'first'));
        trigVal = trigger{i}(trigPos);

        validCount = validCount + 1;
        tmpEpochs(:,validCount) = rec(s:e);
        tmpSigns(validCount,1) = sign(trigVal);
    end

    epochs{i} = tmpEpochs(:,1:validCount);
    epochSigns{i} = tmpSigns(1:validCount);

    fprintf('Type %s: kept %d active epochs, +1=%d, -1=%d\n', ...
        types{i}, validCount, ...
        sum(epochSigns{i}>0), sum(epochSigns{i}<0));
end






avgEpoch = cell(1,4);
alignedEpochs = cell(1,4);
epochsize=info.epochSize
alignStart = round(0.1 * epochSize);
alignEnd   = round(0.45 * epochSize);
maxLag     = round(0.03 * epochSize);

for i = 1:4
    X = epochs{i};
    sgn = epochSigns{i};

    XkeepAll = [];

    for pol = [-1 1]
        Xp = X(:, sgn == pol);

        if isempty(Xp) || size(Xp,2) < 3
            continue;
        end

        ref = median(Xp, 2);

        for iter = 1:3
            X_aligned = zeros(size(Xp));
            refWin = ref(alignStart:alignEnd);

            for k = 1:size(Xp,2)
                x = Xp(:,k);
                xWin = x(alignStart:alignEnd);

                [xc, lags] = xcorr(xWin, refWin, maxLag, 'coeff');
                [~, idxMax] = max(abs(xc));
                lag = lags(idxMax);

                if lag > 0
                    X_aligned(:,k) = [x(lag+1:end); zeros(lag,1)];
                elseif lag < 0
                    s = -lag;
                    X_aligned(:,k) = [zeros(s,1); x(1:end-s)];
                else
                    X_aligned(:,k) = x;
                end
            end

            Xp = X_aligned;
            ref = median(Xp, 2);
        end

        % score epochs in this polarity group
        refWin = ref(alignStart:alignEnd);
        scores = zeros(1, size(Xp,2));

        for k = 1:size(Xp,2)
            x = Xp(:,k);
            xWin = x(alignStart:alignEnd);

            [xc, ~] = xcorr(xWin, refWin, maxLag, 'coeff');
            scores(k) = max(abs(xc));
        end

        keep = scores >= prctile(scores, 60);   % keep top 40%
        Xkeep = Xp(:, keep);

        % bring negative group to common polarity before mixing
        if pol == -1
            Xkeep = -Xkeep;
        end

        XkeepAll = [XkeepAll, Xkeep]; %#ok<AGROW>
    end

    alignedEpochs{i} = XkeepAll;
    avgEpoch{i} = median(XkeepAll, 2);

    fprintf('Type %s: final kept epochs = %d\n', types{i}, size(XkeepAll,2));
end
% figure;
% histogram(epochScores{2}, 20);
% grid on;
% title('Epoch quality scores - Type A');
% xlabel('Score');
% ylabel('Count');
% 
% figure;
% subplot(2,1,1)
% plot(epochs{2}(:,1:min(20,size(epochs{1},2))));
% title('Original epochs - Type A');
% grid on;
% 
% subplot(2,1,2)
% plot(alignedEpochs{2}(:,1:min(20,size(alignedEpochs{1},2))));
% title('Selected aligned epochs - Type A');
% grid on;




% figure;
% subplot(2,1,1)
% plot(epochs{1}(:,1:min(20,size(epochs{1},2))));
% title('Before alignment - first epochs of type A');
% grid on;
% 
% subplot(2,1,2)
% plot(alignedEpochs{1}(:,1:min(20,size(alignedEpochs{1},2))));
% title('After alignment - first epochs of type A');
% grid on;








% -------- align the four averaged stimulus responses to each other --------
avgAligned = avgEpoch;

ref = avgEpoch{1};          % use A as reference
crossMaxLag = 20;           % try 10, 20, 30

for i = 2:4
    x = avgEpoch{i};

    [xc, lags] = xcorr(x(alignStart:alignEnd), ref(alignStart:alignEnd), ...
                       crossMaxLag, 'coeff');
    [~, idxMax] = max(abs(xc));
    lag = lags(idxMax);

    if lag > 0
        avgAligned{i} = [x(lag+1:end); zeros(lag,1)];
    elseif lag < 0
        s = -lag;
        avgAligned{i} = [zeros(s,1); x(1:end-s)];
    else
        avgAligned{i} = x;
    end

    fprintf('Cross-align type %s to A: lag = %d\n', types{i}, lag);
end

% -------- now build OAE from cross-aligned A/B/C/D averages --------
oae_est = avgAligned{1} + avgAligned{2} + avgAligned{3} - 3*avgAligned{4};
oae_est = oae_est(:);
oae_est = oae_est - mean(oae_est); 
oae_est = bandpass(oae_est, [800 4500],fs);









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