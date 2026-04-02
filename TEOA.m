clear; clc;

patient = '194134';
[rec, fs] = audioread(['Patient_' patient '_rec.wav']);
trigger = cell(1,4);
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
    triggerIdx{i} = find(trigger{i} > 0);   
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