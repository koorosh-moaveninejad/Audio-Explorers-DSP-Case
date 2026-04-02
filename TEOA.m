clear; clc;

patient = '194151';
basePath = '/Users/kourosh/Desktop/University/Self Study/Audio-Explorers20206/Diagnostics DSP/Patient Data';
patientFolder = fullfile(basePath, ['patient_' patient]);

[rec, fs] = audioread(fullfile(patientFolder, ['Patient_' patient '_rec.wav']));

active  = cell(1,4);
trigger = cell(1,4);
types = {'A','B','C','D'};

for i = 1:4
    trigger{i} = audioread(fullfile(patientFolder, ['Patient_' patient '_trigger' types{i} '.wav']));
    active{i}  = audioread(fullfile(patientFolder, ['Patient_' patient '_active'  types{i} '.wav']));
end

info = jsondecode(fileread(fullfile(patientFolder, ['Patient_' patient '_info.json'])));
epochSize = info.epochSize;


triggerIdx = cell(1,4);

for i = 1:4
    trigMask = trigger{i} ~= 0;
    actMask  = active{i} > 0.5;

    trigStarts = find(diff([0; trigMask]) == 1);

    actStarts = find(diff([0; actMask]) == 1);
    actEnds   = find(diff([actMask; 0]) == -1);

    validStarts = [];

    for k = 1:length(trigStarts)
        s = trigStarts(k);
        inside = any(s >= actStarts & s <= actEnds);

        if inside
            validStarts(end+1,1) = s; %#ok<AGROW>
        end
    end

    triggerIdx{i} = validStarts;
end

disp(length(triggerIdx{1}))
disp(length(triggerIdx{2}))
disp(length(triggerIdx{3}))
disp(length(triggerIdx{4}))


% 
% are the triggers equally spaced?
% does each trigger define the start of one 882-sample epoch?
% do active masks agree with that?








epochs = cell(1,4);
for i = 1:4
    idx = triggerIdx{i};
    n = length(idx);
    tmp = zeros(epochSize, n);
    validCount = 0;
    for k = 1:n
        s = idx(k);
        e = s + epochSize - 1;
        if e <= length(rec)
            validCount = validCount + 1;
            tmp(:, validCount) = rec(s:e);
        end
    end

    epochs{i} = tmp(:,1:validCount);
end


avgEpoch = cell(1,4);
alignedEpochs = cell(1,4);

maxLag = 20;   % you can tune this later

for i = 1:4
    X = epochs{i};   % size: epochSize x number_of_epochs
    ref = median(X,2);   % much better reference

    X_aligned = zeros(size(X));

    for k = 1:size(X,2)
        x = X(:,k);

        [xc, lags] = xcorr(x, ref, maxLag, 'coeff');
        [~, idxMax] = max(xc);
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

    alignedEpochs{i} = X_aligned;
    avgEpoch{i} = median(X_aligned, 2);
end
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





oae_est = avgEpoch{1} + avgEpoch{2} + avgEpoch{3} - 3*avgEpoch{4};
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
temp = templates.(names{4});
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



% templates = jsondecode(fileread('lostOaes.json'));
% names = fieldnames(templates);
% 
% oae_est = oae_est(:);
% 
% if norm(oae_est) == 0
%     error('oae_est is all zeros');
% end
% 
% oae_n = oae_est / norm(oae_est);
% 
% bestCorr = -inf;
% bestName = '';
% bestTemplate = [];
% 
% for k = 1:length(names)
%     temp = templates.(names{k});
%     temp = temp(:);
% 
%     % resample template to same length as estimated OAE
%     temp_rs = resample(temp, length(oae_est), length(temp));
% 
%     % normalize
%     temp_n = temp_rs / norm(temp_rs);
% 
%     % compare
%     corrVal = dot(oae_n, temp_n);
% 
%     fprintf('%s -> corr = %.3f\n', names{k}, corrVal);
% 
%     if corrVal > bestCorr
%         bestCorr = corrVal;
%         bestName = names{k};
%         bestTemplate = temp_n;
%     end
% end
% 
% fprintf('\nBest match: %s (corr = %.3f)\n', bestName, bestCorr);
% 
% figure;
% plot(oae_n, 'LineWidth', 1.5);
% hold on;
% plot(bestTemplate, 'LineWidth', 1.5);
% grid on;
% legend('Estimated OAE', ['Best template: ' bestName]);
% title(sprintf('Best match: %s (corr = %.3f)', bestName, bestCorr));
% xlabel('Sample');
% ylabel('Normalized amplitude');