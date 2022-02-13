%% Calculates scattering for a single device - in terms of powers and efficiencies
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
clear; clc; close all;

%% Sweep Parameters
thetaOrig = 0;
thetaVals = [thetaOrig-0:1.25:thetaOrig+0];
% thetaVals = [thetaOrig-15:0.625:thetaOrig+15];
peakInd = find(thetaVals==thetaOrig);

%% Generate Corresponding Name of File
sourceConfig = 'gauss';
dp = 3;
file{1} = [pwd,'\'];%,num2str(thetaOrig,'%.1f'),'_inverse_design\'];
file{2} = 'sortspecdata';
file{4} = ['optang',num2str(thetaOrig,'%.1f'),'_th'];
file{3} = [sourceConfig];%, '_xpol'];

fn = formFileName(file,thetaOrig,1);%dp);

%% Main Data-Gathering For Loop
for k = [1:length(thetaVals)]
    theta = thetaVals(k);
    fn = formFileName(file,theta,dp);
    try
        load([fn,'.mat']);
    catch ME
        load([fn(1:end-length(file{4})-1),'.mat']);
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
    
    
    %% Set up Figure Plotting
    sfn = {'device_RTA','exit_scattering'};       % Savefile Names
    realBool = 1;
    
    %% Plot RTA
    fig = figure; hold on;
    xlabel('Wavelength (um)');
    ylabel('Normalised Power');
    
    ov = zeros(15,1);
    if realBool
        sfn{1} = [sfn{1},'_real'];
    else
        ov = [0 0 0.06];
    end

    intensity = 230;
    plot(wlVals, dev0.ref./dev0.P_in, 'o-', 'Color',[0 0 intensity]./255, 'DisplayName', 'Reflection');
    plot(wlVals, (dev0.P_out)./dev0.P_in, ...
        'o-', 'Color',[0 intensity 0]./255, 'DisplayName', 'Transmission_{bottom}');
    plot(wlVals, (dev0.sm0+dev0.sm1+dev0.sm2+dev0.sm3)./dev0.P_in, ...
        'o-', 'Color',[0 intensity/2 0]./255, 'DisplayName', 'Transmission_{side}');
    plot(wlVals, (dev0.abs)./dev0.P_in + ov(3), 'o-', 'Color',[intensity 0 0]./255, 'DisplayName', 'Absorption');
%     plot(wlVals, (dev0.ref+dev0.abs+dev0.P_out+dev0.sm0+dev0.sm1+dev0.sm2+dev0.sm3)./dev0.P_in, ...
%         'o-', 'Color',[0 0 0]./255, 'DisplayName', 'Sum');

    if ~realBool
        ylim([0,1]);
    end
    leg = legend('Location', 'east');
    title(['RTA of Device: normalized against Input Aperture']);
    
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
    
    %% Plot Scattering Efficiencies
    fig = figure; hold on;
    xlabel('Wavelength (um)');
    ylabel('Normalised Power');
    
    ov = zeros(15,1);
    if realBool
        sfn{2} = [sfn{2},'_real'];
    end

    intensity = 230;
    plot(wlVals, dev0.P_fp./dev0.P_out, 'o-', 'Color',[0 intensity 0]./255, 'DisplayName', 'Focal Region');
    plot(wlVals, dev0.sm0./dev0.P_out, 'o-', 'Color',[intensity/4 0 0]./255, 'DisplayName', 'Side N');
    plot(wlVals, dev0.sm1./dev0.P_out, 'o-', 'Color',[2*intensity/4 0 0]./255, 'DisplayName', 'Side S');
    plot(wlVals, dev0.sm2./dev0.P_out, 'o-', 'Color',[3*intensity/4 0 0]./255, 'DisplayName', 'Side E');
    plot(wlVals, dev0.sm3./dev0.P_out, 'o-', 'Color',[intensity 0 0]./255, 'DisplayName', 'Side W');
    plot(wlVals, (dev0.P_out - dev0.P_fp)./dev0.P_out, ...
        'o-', 'Color',[0 0 intensity]./255, 'DisplayName', 'Oblique Scattering');
    %plot(wlVals, ((dev0.P_spill - dev0.P_miss - dev0.P_fp)+dev0.P_fp+dev0.sm0+dev0.sm1+dev0.sm2+dev0.sm3)./dev0.P_out);

    if ~realBool
        ylim([0,inf]);
    end
    leg = legend('Location', 'east');
    title(['Scattering at Focal Plane: normalized against Exit Aperture']);
    
    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
      lines(i).LineWidth = 2.0;
      lines(i).MarkerSize = 2.0;
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    %set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
    set(gcf,'position',[0 0 1920 1440]);
    
    exportgraphics(gca,[sfn{2},'.png']);
    saveas(gca,[sfn{2},'.fig']);
    %close all;

end

%% Functions
function fn = formFileName(file,theta,dp)
    switch dp
        case 0
            thetaStr = num2str(theta,'%.0f');
        case 1
            thetaStr = num2str(theta,'%.1f');
        case 2
            thetaStr = num2str(theta,'%.2f');
        case 3
            thetaStr = num2str(theta,'%.3f');
        otherwise
            thetaStr = num2str(theta,'%d');
    end

    fn = [];
    for i = 1:length(file)
        if file{i}(end-2:end) == '_th'
           file{i} = [file{i},thetaStr]; 
        end
        
        if i <= 2
            fn = [fn,file{i}];
        else
            fn = [fn,'_',file{i}];
        end
    end
end

function [power_lambda] = integrateOverSpace(monitor_name,variable,theta)
    Emag = magnitudeField(monitor_name,variable);
    
    if variable=='E'
        Emag_lambda = squeeze(sum(Emag,[1,2,3]));
        intensity_lambda = 0.5*3e8*8.85418782e-12*1*cos(theta)*Emag_lambda.^2;
        power_lambda = intensity_lambda*(9e-5*9e-5);
    else
        power_lambda = abs(Emag);
    end
end

function [mag_xyz_lambda] = magnitudeField(monitor_name,variable)
    monitor_name = process_dataset(monitor_name,variable);

    field = monitor_name.(variable);
    if variable=='E' || variable=='H'
        fieldX = field(:,:,:,:,1);
        fieldY = field(:,:,:,:,2);
        fieldZ = field(:,:,:,:,3);
        mag_xyz_lambda = sqrt(abs(fieldX).^2+abs(fieldY).^2+abs(fieldZ).^2);
    else
        mag_xyz_lambda = field;
    end
end

function [E] = process_dataset(monitor,field)
    monitor.vars.(field) = monitor.data;
    E = monitor.vars;
end

function [E_spatial] = reshape_spatial(monitor_name,variable,field_component,wlVal)
%     Deprecated because even though this is probably right I'd rather not
%     deal with this: Ian Foo, 20210924
%     component: x,y,z = 1,2,3
%     wlVal = index in the wavelength value array
    s1 = length(monitor_name.x);
    s2 = length(monitor_name.y);
    s3 = length(monitor_name.z);
    
    E = monitor_name.(variable)(:,field_component,wlVal);
    
    E_spatial = reshape(E,[s1 s2 s3]);
end