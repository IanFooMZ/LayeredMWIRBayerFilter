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

sweepVals = [sweepOrig-0:2:sweepOrig+0];
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
    dev0.lambda = wlVals;
	dev0.x = out.E_fp0.vars.x;
	dev0.y = out.E_fp0.vars.y;
   
    % Power hitting Focal Region
    dev0.E_fp = abs(out.E_fp0.vars.E);
    % Power in entire FDTD focal plane
    dev0.E_spill = abs(out.E_sp0.vars.E);
    % At this point have not added aperture to attenuate power missing the
    % device yet
    
    % Power hitting each focal monitor
    dev0.E_tm0 = abs(out.E_tm0.vars.E);
    dev0.E_tm1 = abs(out.E_tm1.vars.E);
    dev0.E_tm2 = abs(out.E_tm2.vars.E);
    dev0.E_tm3 = abs(out.E_tm3.vars.E);

    devCell{m} = dev0;

	%% Set up Figure Plotting
    realBool = 0;
    wlPlotVals = [3.5 4.5 5.5]; % in um
    
    fig = monitorPattern(devCell{m}, 'E_spill', wlPlotVals, 'norm', realBool);
end

function fig = monitorPattern(struct,field, wlValues,vectorType, realBool)

    %% Set up Figure Plotting
    realStr = '';
    [d,idx] = min(abs(struct.lambda-wlValues));
    for index = idx
        plotWL = struct.lambda(index);
        
        % Figure Properties
        fig = figure; ax = gca; hold on;
        xlabel('y'); ylabel('x');
        title({['Scatt. Pattern in Focal Plane'],...
            ['\lambda: ',num2str(plotWL,'%.1f'),' um']});
        colormap jet; colorbar;
        set(gcf,'position',[322 36 600 600]);
        
        
        if realBool
            realStr = '_real';
        end
        sfn = {[field, '_', vectorType, '_',...
            num2str(plotWL,'%.1f'),'um',...
            realStr]};		% Savefile Names
        
        %% Process field and rearrange array
        plotData = struct.(field);
        % Get vector dimension
        switch vectorType
            case 'x'
                plotData = plotData(:,:,:,:,1);
            case 'y'
                plotData = plotData(:,:,:,:,2);
            case 'z'
                plotData = plotData(:,:,:,:,3);
            otherwise
                fieldX = plotData(:,:,:,:,1);
                fieldY = plotData(:,:,:,:,2);
                fieldZ = plotData(:,:,:,:,3);
                plotData = sqrt(abs(fieldX).^2+abs(fieldY).^2+abs(fieldZ).^2);
        end
        % Take the slice at the wavelength of the interest
        plotData = plotData(:,:,:,index);
        % Squeeze z dimension
        plotData = squeeze(plotData);
        
        %% Plot Scattering Pattern at the Wavelength
        szData = size(plotData);
        [X,Y] = meshgrid(struct.y,struct.x);
        h = surf(X,Y,plotData);
        
        
        % Figure Post-Processing
        set(h,'LineStyle','none');
        view(90,-90); 
        axis square;
        exportgraphics(gca,[sfn{1},'.png']);
        %saveas(gca,[sfn{1},'.fig']);
        close all;
    
    end
end
