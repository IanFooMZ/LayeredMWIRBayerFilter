%% Calculates scattering for a range of devices - in terms of powers and efficiencies
% Definition: 
% Note 59 in OneNote PhD Notebook: https://bit.ly/3BeJXoE
% Note 60 in OneNote PhD Notebook: https://bit.ly/3ozgpOP
% 	1. Source power
% 	2. Transmission Monitor -> Reflection (normalized against 1)
% 	3. Incident Plane Monitor -> How much doesn't even pass through the device at all (n.a. 2)
% 	4. Side scattering (n.a. 3)
% 	5. Exit Aperture (n.a. 3)
% 	6. Focal Region (n.a. 5 OR 6+7)
% 		a. Has to be n.a. 6+7 because we haven't attenuated the power that skips the device yet
% 	7. Outside Focal Region (n.a. 5 OR 6+7)
% 		a. Has to be n.a. 6+7 b/c above
%         b. See if it increases when the FDTD region does

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

wlPickInd = 1; %4.5um - we have to slice the P_arrs at some wavelength, go for green
% blue is 1; red is end

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

save('temp');
% load('temp.mat');

%% Collect specific statistics for the entire array of devices in their respective
%% 2D arrays: (sweepVal, wavelength) are the dimensions
Parr_ip = fc.extractParamVals(devCell,'P_in');
Parr_fp = fc.extractParamVals(devCell,'P_fp');
Parr_spill = fc.extractParamVals(devCell,'P_spill');
Parr_tm0 = fc.extractParamVals(devCell,'P_tm0');
Parr_tm1 = fc.extractParamVals(devCell,'P_tm1');
Parr_tm2 = fc.extractParamVals(devCell,'P_tm2');
Parr_tm3 = fc.extractParamVals(devCell,'P_tm3');

Parr_ref = fc.extractParamVals(devCell,'ref');
Parr_out = fc.extractParamVals(devCell,'P_out');
Parr_sm0 = fc.extractParamVals(devCell,'sm0');
Parr_sm1 = fc.extractParamVals(devCell,'sm1');
Parr_sm2 = fc.extractParamVals(devCell,'sm2');
Parr_sm3 = fc.extractParamVals(devCell,'sm3');
Parr_abs = fc.extractParamVals(devCell,'abs');

%% Take slice along wavelength dimension (designed peak),
Parr_fp = Parr_fp(:,wlPickInd);
Parr_ip = Parr_ip(:,wlPickInd);
Parr_ref = Parr_ref(:,wlPickInd);
Parr_out = Parr_out(:,wlPickInd);
Parr_sm0 = Parr_sm0(:,wlPickInd);
Parr_sm1 = Parr_sm1(:,wlPickInd);
Parr_sm2 = Parr_sm2(:,wlPickInd);
Parr_sm3 = Parr_sm3(:,wlPickInd);
Parr_abs = Parr_abs(:,wlPickInd);



%% Set up Figure Plotting
realBool = 0; realStr = ''; ov = zeros(15,1);

% Figure Properties
fig = figure; ax = gca; hold on;
xlabel('Sidewall Thickness (um)');
ylabel('Normalized Power');
%xlim([sweepVals(4) sweepVals(end-3)]); 
%ylim([0.4,1.05]);
legend = legend('Location', 'northeast');
title(['Sidewall Thickness Sweep: Optimized at ',num2str(thetaOrig,"%.1f"),'°']);

if realBool
    realStr = '_real';
else
    %ov = fc.findOffset(peakInd,{Parr_tm0,Parr_tm1,Parr_tm2});
end
sfnStr = [fc.compileParamStr(paramCell(1:end-3)),...
    realStr];
sfn = {['device_RTA', sfnStr],...
    ['exit_scattering',sfnStr]};		% Savefile Names

% xline(thetaOrig,'-','HandleVisibility','off','Color','#505050' ...
%     ,'LineWidth', 2.0);
% yline(0.5,'--','HandleVisibility','off','Color','#505050' ...
%     ,'LineWidth', 2.0);

%% Plot Sorting Spectrum
% intensity = 230;
% plot(sweepVals+ov(1)*dswp,Parr_tm0,'o-','Color',[0 0 intensity]./255,'DisplayName','Blue');
% plot(sweepVals+ov(2)*dswp,Parr_tm1,'o-','Color',[0 intensity 0]./255,'DisplayName','Green, x-pol');
% plot(sweepVals+ov(3)*dswp,Parr_tm2,'o-','Color',[intensity 0 0]./255,'DisplayName','Red');
% % plot(sweepVals+ov(4)*dswp,Parr_tm3,'Color',1/255*[40,94,25],'DisplayName','Green,y-pol');

intensity = 230;
% YOU GOT TO TAKE THE THING AT THE WAVELENGTH OF INTEREST
%plot(sweepVals, Parr_ref./Parr_ip, 'o-', 'Color',[0 0 intensity]./255, 'DisplayName', 'Reflection');
plot(sweepVals, Parr_out./Parr_ip, ...
    'o-', 'Color',[0 intensity 0]./255, 'DisplayName', 'Transmission_{bottom}');
plot(sweepVals, (Parr_sm0+Parr_sm1+Parr_sm2+Parr_sm3)./Parr_ip, ...
    'o-', 'Color',[0 intensity/2 0]./255, 'DisplayName', 'Transmission_{side}');
%plot(sweepVals, Parr_abs./Parr_ip, 'o-', 'Color',[intensity 0 0]./255, 'DisplayName', 'Absorption');
%     plot(wlVals, (dev0.ref+dev0.abs+dev0.P_out+dev0.sm0+dev0.sm1+dev0.sm2+dev0.sm3)./dev0.P_in, ...
%         'o-', 'Color',[0 0 0]./255, 'DisplayName', 'Sum');

% Figure Post-Processing
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
  lines(i).MarkerSize = 2.0;
end
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
set(gcf,'position',[0 0 1920 1440]);

transBtm = Parr_out./Parr_ip;
transSide = (Parr_sm0+Parr_sm1+Parr_sm2+Parr_sm3)./Parr_ip;
save('testSide','transBtm','transSide');

exportgraphics(gca,[sfn{1},'.png']);
saveas(gca,[sfn{1},'.fig']);
%close all;