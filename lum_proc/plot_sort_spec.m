%% Calculates sorting efficiency for each quadrant
% Definition: Fraction of incident power reaching target quadrant
% References: 
% https://support.lumerical.com/hc/en-us/articles/360034409554
% https://kx.lumerical.com/t/transforming-datasets-as-structures-to-matlab/2576/5
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

%% Generate Corresponding Name of File
compileFN;  %addpath(pwd);
%% Store data for each run in a cell array
devCell = cell(length(sweepVals),1);

%% Main Data-Gathering For Loop
for m = [1:length(sweepVals)]
    swp = sweepVals(m); paramCell{end-1} = swp;
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
%     dev0.sm0 = abs(out.E_sm0.power);
%     dev0.sm1 = abs(out.E_sm1.power);
%     dev0.sm2 = abs(out.E_sm2.power);
%     dev0.sm3 = abs(out.E_sm3.power);
    % Exit Aperture
    dev0.P_out = abs(out.E_eam.power);
    % Power absorbed by Device
    %dev0.abs = (dev0.P_in - dev0.ref - dev0.sm0 - dev0.sm1 - dev0.sm2 - dev0.sm3 - dev0.P_out);
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

	%% Set up Figure Plotting
    transPlot = 1; transStr = '';
    realBool = 0; realStr = ''; ov = zeros(15,1);
    
    % Figure Properties
    fig = figure; ax = gca; hold on;
	xlabel('Wavelength (um)');
    ylabel('Sorting Efficiency');
    ylim([0,0.9]);
    leg = legend('Location', 'north');
    title({['Spectrum in Each Quadrant'],...
        ['Incident Angle ',num2str(swp), ...
        '°, Optimized for ',num2str(thetaOrig),'°']});
    
	if transPlot
		transStr = '_trans';
        ylim([0.1,0.8]);
        leg = legend('Location', 'northwest');
    end
    if realBool
        realStr = '_real';
    else
        ov = [0 0 0.06];
    end
    sfn = {['sortspec', transStr, ...
        erase(fn(length([pwd,'\'])+1:end),file{2}),...
        realStr]};		% Savefile Names

    %% Plot Sorting Spectrum
    intensity = 230;
    plot(wlVals,dev0.P_tm0./dev0.P_fp,'o-','Color',[0 0 intensity]./255,'DisplayName','Blue');
    plot(wlVals,dev0.P_tm1./dev0.P_fp,'o-','Color',[0 intensity 0]./255,'DisplayName','Green, x-pol');
    plot(wlVals,dev0.P_tm2./dev0.P_fp,'o-','Color',[intensity 0 0]./255,'DisplayName','Red');
    plot(wlVals,dev0.P_tm3./dev0.P_fp,'Color',1/255*[40,94,25],'DisplayName','Green, y-pol');
    if transPlot
        yyaxis right; ylabel('Transmission'); ax.YAxis(2).Color = [51 0 0]./255;
        plot(wlVals,dev0.P_fp./dev0.P_in,'o-','Color',[51 0 0]./255,'DisplayName','Transmission');
    end

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
	%saveas(gca,[sfn{1},'.fig']);
    close all;

end