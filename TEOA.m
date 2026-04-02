clear; clc;

patient = '194008';
basePath = '/Users/kourosh/Desktop/University/Self Study/Audio-Explorers20206/Diagnostics DSP/Patient Data';
filePath = fullfile(basePath, ['patient_' patient], ['Patient_' patient '_rec.wav']);
[rec, fs] = audioread(filePath);


active  = cell(1,4);
types = {'A','B','C','D'};
for i = 1:4
    trigger{i} = audioread(['Patient_' patient '_trigger' types{i} '.wav']);
    active{i}  = audioread(['Patient_' patient '_active'  types{i} '.wav']);
end
info = jsondecode(fileread(['Patient_' patient '_info.json']));
epochSize = info.epochSize;


triggerIdx = cell(1,4);

for i = 1:4
    trigMask = trigger{i} > 0;
    actMask  = active{i} > 0.5;

    % trigger starts: rising edges of trigger
    trigStarts = find(diff([0; trigMask]) == 1);

    % active regions: start/end of active blocks
    actStarts = find(diff([0; actMask]) == 1);
    actEnds   = find(diff([actMask; 0]) == -1);

    validStarts = [];

    for k = 1:length(trigStarts)
        s = trigStarts(k);

        % keep trigger only if it falls inside an active block
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

for i = 1:4
    avgEpoch{i} = mean(epochs{i}, 2);
end

figure;
for i = 1:4
    subplot(4,1,i);
    plot(avgEpoch{i});
    title(['Average epoch ' types{i}]);
end











oae_est = avgEpoch{1} + avgEpoch{2} + avgEpoch{3} - 3*avgEpoch{4};
t = (0:epochSize-1)/fs*1000; % ms

figure;
plot(t, oae_est);
xlabel('Time (ms)');
ylabel('Amplitude');
title(['Estimated OAE for patient ' patient]);
grid on;

templates = jsondecode(fileread('lostOaes.json'));
names = fieldnames(templates);

oae_est = oae_est(:);

if norm(oae_est) == 0
    error('oae_est is all zeros');
end

oae_n = oae_est / norm(oae_est);

bestCorr = -inf;
bestName = '';
bestTemplate = [];

for k = 1:length(names)
    temp = templates.(names{k});
    temp = temp(:);

    % resample template to same length as estimated OAE
    temp_rs = resample(temp, length(oae_est), length(temp));

    % normalize
    temp_n = temp_rs / norm(temp_rs);

    % compare
    corrVal = dot(oae_n, temp_n);

    fprintf('%s -> corr = %.3f\n', names{k}, corrVal);

    if corrVal > bestCorr
        bestCorr = corrVal;
        bestName = names{k};
        bestTemplate = temp_n;
    end
end

fprintf('\nBest match: %s (corr = %.3f)\n', bestName, bestCorr);

figure;
plot(oae_n, 'LineWidth', 1.5);
hold on;
plot(bestTemplate, 'LineWidth', 1.5);
grid on;
legend('Estimated OAE', ['Best template: ' bestName]);
title(sprintf('Best match: %s (corr = %.3f)', bestName, bestCorr));
xlabel('Sample');
ylabel('Normalized amplitude');