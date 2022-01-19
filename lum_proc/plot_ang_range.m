%% Calculates sorting efficiency spectrum for each quadrant
% Definition: Fraction of incident power reaching target quadrant
% References: 
% https://support.lumerical.com/hc/en-us/articles/360034409554
% https://kx.lumerical.com/t/transforming-datasets-as-structures-to-matlab/2576/5

%% And then plots it according to the peak sorting efficiency at each angle
clear; clc; close all;
fc = functionsContainer;

%% Sweep Parameters
%% Only uncomment this block if parameters differ from those in the compileFN.m
% paramCell = {'th',5,'%.1f',...
%             'bRad',5,'%.3f'};
% Format is (variable name, variable start value, format string).
% Last one is the sweep variable
% sweepOrig = paramCell{end-1};
% sweepVals = [sweepOrig-0:2:sweepOrig+40];

%% Generate Corresponding Name of Files
compileFN;  %addpath(pwd);

idcs = strfind(pwd,'\'); mydir = pwd; newdir = mydir(1:idcs(end)-1);
compileFileName = [newdir,'\','angRange_superimp',sourceConfig,'_',paramCell{end-2},'.mat'];
saveVarName = 'th5bRad_test';
% saveVarName = fc.compileParamStr(paramCell(1:end-3));
% saveVarName(regexp(saveVarName,'[.,0,_]'))=[];

%% Store data for each run in a cell array
devCell = cell(length(sweepVals),1);

%% Main Data-Gathering For Loop
for m = [1:length(sweepVals)]
    swp = sweepVals(m); paramCell{end-1} = swp;
    disp(['(',num2str(m),'/',num2str(length(sweepVals)),')   ',...
        paramCell{end-2},': ',num2str(swp)]);
    fn = fc.formFileName(file,paramCell);
    try
        load([fn,'.mat']);
    end
    wlVals = 1e6*(out.lambda);
	
	%% Source Power Spectrum
    sp = out.sourcepower;
    %% Central Device Specs (all following measurements are in power)
    dev0 = struct;
    dev0.sp = sp;
    % Reflected from Device
    dev0.ref = (sp-abs(out.E_stm0.power));
    % Incident Plane Monitor
    dev0.P_in = abs(out.E_iam.power);
    % Power that Misses Device
    dev0.P_miss = (dev0.sp - dev0.ref) - dev0.P_in;
    % Side Scattering
    dev0.sm0 = abs(out.E_sm0.power);
    dev0.sm1 = abs(out.E_sm1.power);
    dev0.sm2 = abs(out.E_sm2.power);
    dev0.sm3 = abs(out.E_sm3.power);
    % Exit Aperture
    dev0.P_out = abs(out.E_eam.power);
    % Power absorbed by Device
    dev0.abs = (dev0.P_in - dev0.ref - dev0.sm0 - dev0.sm1 - dev0.sm2 - dev0.sm3 - dev0.P_out);
    % Power hitting Focal Region
    dev0.P_fp = abs(out.E_fp0.power);
    % Power in entire FDTD focal plane
    dev0.P_spill = abs(out.E_sp0.power);
    % At this point have not added aperture to attenuate power missing the
    % device yet
    
    % Power hitting each focal monitor
    dev0.P_tm0 = abs(out.E_tm0.power);
    dev0.P_tm1 = abs(out.E_tm1.power);
    dev0.P_tm2 = abs(out.E_tm2.power);
    dev0.P_tm3 = abs(out.E_tm3.power);

    devCell{m} = dev0;
end

%% Collect specific statistics for the entire array of devices in their respective
%% 2D arrays: (sweepVal, wavelength) are the dimensions
Parr_ip = fc.extractParamVals(devCell,'P_in');
Parr_fp = fc.extractParamVals(devCell,'P_fp');
Parr_spill = fc.extractParamVals(devCell,'P_spill');
Parr_tm0 = fc.extractParamVals(devCell,'P_tm0');
Parr_tm1 = fc.extractParamVals(devCell,'P_tm1');
Parr_tm2 = fc.extractParamVals(devCell,'P_tm2');
Parr_tm3 = fc.extractParamVals(devCell,'P_tm3');


%% Calculate contrast of each bin: defined as the highest intensity bin /
%% 2nd highest intensity bin
contrast_tm0 = fc.calcContrast(Parr_tm0,Parr_tm1,Parr_tm2,Parr_tm3,peakInd);
contrast_tm1 = fc.calcContrast(Parr_tm1,Parr_tm0,Parr_tm2,Parr_tm3,peakInd);
contrast_tm2 = fc.calcContrast(Parr_tm2,Parr_tm1,Parr_tm0,Parr_tm3,peakInd);
contrast_tm3 = fc.calcContrast(Parr_tm3,Parr_tm1,Parr_tm2,Parr_tm0,peakInd);


%% Normalize along focal plane,
%% Take Maximum along wavelength dimension (designed peak),
%% Normalize to max. value
Parr_spill = fc.findMax((Parr_spill-Parr_fp)./(Parr_spill), peakInd);
Parr_tm0 = fc.findMax(Parr_tm0./Parr_fp, peakInd);
Parr_tm1 = fc.findMax(Parr_tm1./Parr_fp, peakInd);
Parr_tm2 = fc.findMax(Parr_tm2./Parr_fp, peakInd);
Parr_tm3 = fc.findMax(Parr_tm3./Parr_fp, peakInd);

%% Save to Compile File as a Struct
% Save to struct
svData = struct;
svData.wlVals = wlVals;
svData.peakInd = peakInd;
svData.sweepVals = sweepVals;
svData.paramCell = paramCell;
svData.Parr_ip = Parr_ip;
svData.Parr_fp = Parr_fp;
svData.Parr_spill = Parr_spill;
svData.Parr_tm0 = Parr_tm0;
svData.Parr_tm1 = Parr_tm1;
svData.Parr_tm2 = Parr_tm2;
svData.Parr_tm3 = Parr_tm3;
svData.ov = fc.findOffset(peakInd,{Parr_tm0,Parr_tm1,Parr_tm2,Parr_tm3});

S.(saveVarName) = svData;
if isfile(compileFileName)
	save(compileFileName,'-struct','S','-append');
else    save(compileFileName,'-struct','S');
end


%% Normalization and Adjustments of Array Values
Parr_tm0 = Parr_tm0./Parr_tm0(peakInd);
Parr_tm1 = Parr_tm1./Parr_tm1(peakInd);
Parr_tm2 = Parr_tm2./Parr_tm2(peakInd);
Parr_tm3 = Parr_tm3./Parr_tm3(peakInd);

dswp = sweepVals(2)-sweepVals(1);

Parr_spill = Parr_spill./max(Parr_spill);
Parr_tm0 = Parr_tm0./max(Parr_tm0);
Parr_tm1 = Parr_tm1./max(Parr_tm1);
Parr_tm2 = Parr_tm2./max(Parr_tm2);


%% Set up Figure Plotting
realBool = 0; realStr = ''; ov = zeros(15,1);

% Figure Properties
fig = figure; ax = gca; hold on;
xlabel('Angle of Incidence (°)');
ylabel('Sorting Efficiency (Normalized)');
xlim([sweepVals(4) sweepVals(end-3)]); 
ylim([0.4,1.05]);
legend = legend('Location', 'northeast');
title(['Angular Range: Optimized at ',num2str(thetaOrig,"%.1f"),'°']);

if realBool
    realStr = '_real';
else
    ov = fc.findOffset(peakInd,{Parr_tm0,Parr_tm1,Parr_tm2});
end
sfn = {['angrange', ...
    fc.compileParamStr(paramCell(1:end-3)),...
    realStr]};		% Savefile Names

xline(thetaOrig,'-','HandleVisibility','off','Color','#505050' ...
    ,'LineWidth', 2.0);
yline(0.5,'--','HandleVisibility','off','Color','#505050' ...
    ,'LineWidth', 2.0);

%% Plot Sorting Spectrum
intensity = 230;
plot(sweepVals+ov(1)*dswp,Parr_tm0,'o-','Color',[0 0 intensity]./255,'DisplayName','Blue');
plot(sweepVals+ov(2)*dswp,Parr_tm1,'o-','Color',[0 intensity 0]./255,'DisplayName','Green, x-pol');
plot(sweepVals+ov(3)*dswp,Parr_tm2,'o-','Color',[intensity 0 0]./255,'DisplayName','Red');
% plot(sweepVals+ov(4)*dswp,Parr_tm3,'Color',1/255*[40,94,25],'DisplayName','Green,y-pol');


% Figure Post-Processing
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
  lines(i).MarkerSize = 2.0;
end
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
set(gcf,'position',[0 0 1920 1440]);

exportgraphics(gca,[sfn{1},'.png']);
saveas(gca,[sfn{1},'.fig']);
%close all;