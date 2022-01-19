%% In the same folder, there will be a list of save data files named according to
% parameters and the last param. will be the variable that is being swept.
fc = functionsContainer;

%% Sweep Parameters
paramCell = {'th',5,'%.1f',...
            'bRad',5,'%.3f'};
% Format is (variable name, variable start value, format string).
% Last one is the sweep variable

sweepOrig = paramCell{end-1};
sweepVals = [sweepOrig+0:2:sweepOrig+40];
peakInd = find(sweepVals==sweepOrig);

thetaOrig = paramCell{2};   % Hardcode this part

%% Generate Corresponding Name of File
sourceConfig = 'gauss';

file{1} = [pwd,'\'];
file{2} = ['sortspecdata'];
file{3} = [sourceConfig];
fn = fc.formFileName(file,paramCell);

% k = 5;
% 
% swp = sweepVals(k); paramCell{end-1} = swp;
% fn = formFileName(file,paramCell);
% try
%     load([fn,'.mat']);
% catch ME
%     load([fn(1:end-length(file{4})-1),'.mat']);
% end
