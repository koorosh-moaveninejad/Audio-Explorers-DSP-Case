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



