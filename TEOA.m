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

