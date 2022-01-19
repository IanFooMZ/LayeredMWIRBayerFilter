clear; clc; close all;

%% Sweep Parameters
%% Only uncomment this block if parameters differ from those in the compileFN.m
% paramCell = {'th',5,'%.1f',...
%             'bRad',5,'%.3f'};
% Format is (variable name, variable start value, format string).
% Last one is the sweep variable
% sweepOrig = paramCell{end-1};
% sweepVals = [sweepOrig-0:2:sweepOrig+40];


%% Generate Corresponding Name of File
compileFN;  %addpath(pwd);

%% Sweep through Each File
for k = [1:length(sweepVals)]
    swp = sweepVals(k); paramCell{end-1} = swp;
    fn = fc.formFileName(file,paramCell);
    disp(['Processing ', paramCell{end-2}, num2str(paramCell{end-1},paramCell{end}),'...']);
    try
        load([fn,'.mat']);
    end

    fieldNms = fieldnames(out);

    for k = 1:length(fieldNms)
        field = fieldNms{k};
        if strcmp(class(out.(field)), 'char')
            rmfield(out,field);
        elseif strcmp(class(out.(field)), 'struct')
            out.(field) = process_dataset(out.(field), field(1));
        end
    end

    save([fn,'.mat'],'out');
end

function [E] = process_dataset(monitor,field)
    if ~any(strcmp(fieldnames(monitor),'data'))
    	E = monitor;
    else
        monitor.vars.(field) = monitor.data;
        monitor = rmfield(monitor,'data');
    	E = monitor;
    end
end
