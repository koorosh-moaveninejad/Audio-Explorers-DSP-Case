clear; clc;

patient = '194143';
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
maxLag = 40;
alignStart = 30;
alignEnd   = min(epochSize, 40);

for i = 1:4
    X = epochs{i};

    ref = median(X, 2);

    for iter = 1:3
        X_aligned = zeros(size(X));

        refWin = ref(alignStart:alignEnd);

        for k = 1:size(X,2)
            x = X(:,k);
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

        X = X_aligned;
        ref = median(X, 2);
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
oae_est = oae_est(:);
oae_est = oae_est - mean(oae_est);
oae_est = bandpass(oae_est, [1000 3900], fs);


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
temp = templates.(names{3});
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






figure;
hold on;

for k = 1:min(10, size(alignedEpochs{1},2))
    x = alignedEpochs{1}(:,k);
    x = x - mean(x);

    N = length(x);
    Y = fft(x);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:floor(N/2))/N;

    plot(f, P1);
end

grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of individual aligned epochs (Type A)');
xlim([0 5000]);

figure;
hold on;


for k = 1:min(10, size(alignedEpochs{2},2))
    x = alignedEpochs{2}(:,k);
    x = x - mean(x);

    N = length(x);
    Y = fft(x);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:floor(N/2))/N;

    plot(f, P1);
end

grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of individual aligned epochs (Type A)');
xlim([0 5000]);