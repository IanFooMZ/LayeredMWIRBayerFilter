%% Calculates sorting efficiency for each quadrant
% Definition: Fraction of incident power reaching target quadrant
% References: 
% https://support.lumerical.com/hc/en-us/articles/360034409554
% https://kx.lumerical.com/t/transforming-datasets-as-structures-to-matlab/2576/5
clear; clc; close all;
thetaOrig = 10;
%thetaVals = [thetaOrig:1.25:thetaOrig+30];
thetaVals = [thetaOrig-15:1.25:thetaOrig+15];
peakInd = find(thetaVals==thetaOrig);
sourceConfig = 'gauss';
dp = 3;

file{1} = [pwd,'\'];%,num2str(thetaOrig,'%.1f'),'_inverse_design\'];
file{2} = 'sortspecdata';
file{4} = ['optang',num2str(thetaOrig,'%.1f'),'_th'];
file{3} = [sourceConfig];%, '_xpol'];

fn = formFileName(file,thetaOrig,1);%dp);

for k = [1:length(thetaVals)]
    theta = thetaVals(k);
%     if theta<0
%         fn = formFileName(file,ceil(theta),dp);
%     else
%         fn = formFileName(file,floor(theta),dp);
%     end
    fn = formFileName(file,theta,dp);
    try
        load([fn,'.mat']);
    catch ME
        load([fn(1:end-length(file{4})-1),'.mat']);
    end
    wlVals = 1e6*E_fm0.lambda;
    
    %% Incident Power at Input Aperture
    Emag_im0 = integrateOverSpace(E_im0,'E',theta);
    Emag_ip0 = integrateOverSpace(E_ip0,'E',theta);
    
    %% Overall Incident Power
    Emag_fp0 = integrateOverSpace(E_fp0,'E',theta);

    %% Incident Power for Each Quadrant
    Emag_tm0 = integrateOverSpace(E_tm0,'E',theta);
    Emag_tm1 = integrateOverSpace(E_tm1,'E',theta);
    Emag_tm2 = integrateOverSpace(E_tm2,'E',theta);
    Emag_tm3 = integrateOverSpace(E_tm3,'E',theta);

    %% E-Field at Each Focal Monitor
%     Emag_fm0 = integrateOverSpace(E_fm0,'E',theta);
%     Emag_fm1 = integrateOverSpace(E_fm1,'E',theta);
%     Emag_fm2 = integrateOverSpace(E_fm2,'E',theta);
%     Emag_fm3 = integrateOverSpace(E_fm3,'E',theta);

    %% Plot Sorting Spectrum
    fig = figure; hold on;
    transPlot = true;
    intensity = 230;
    if transPlot
        plot(wlVals,Emag_fp0./Emag_ip0,'o-','Color',[0 0 0]./255,'DisplayName','Transmission');
    end
    plot(wlVals,Emag_tm0./Emag_fp0,'o-','Color',[0 0 intensity]./255,'DisplayName','Blue');
    plot(wlVals,Emag_tm1./Emag_fp0,'o-','Color',[0 intensity 0]./255,'DisplayName','Green, x-pol');
    plot(wlVals,Emag_tm2./Emag_fp0,'o-','Color',[intensity 0 0]./255,'DisplayName','Red');
    plot(wlVals,Emag_tm3./Emag_fp0,'Color',1/255*[40,94,25],'DisplayName','Green, y-pol');

    xlabel('Wavelength (um)');
    ylabel('Sorting Efficiency');
    if transPlot
        ylim([0.1,0.8]);
        leg = legend('Location', 'northwest');
    else
    	ylim([0.1,0.55]);
        leg = legend('Location', 'north');
    end
    title({['Spectrum in Each Quadrant'],['Incident Angle ',num2str(theta) ...
        '°, Optimized for ',num2str(thetaOrig),'°']});
    
    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
      lines(i).LineWidth = 2.0;
      lines(i).MarkerSize = 2.0;
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    %set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
    set(gcf,'position',[0 0 1920 1440]);
    
    if transPlot
        saveFileName = ['sortspec_',sourceConfig,'_trans_optang', num2str(thetaOrig) ...
        ,'_th',num2str(theta),'.png'];
    else
        saveFileName = ['sortspec_',sourceConfig,'_optang', num2str(thetaOrig) ...
        ,'_th',num2str(theta),'.png'];
    end
    exportgraphics(gca,saveFileName);

    close all;

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

function [Emag_lambda] = integrateOverSpace(monitor_name,variable,theta)
    Emag = magnitudeE(monitor_name,variable);
    Emag_lambda = squeeze(sum(Emag,1));
    
    power_lambda = 0.5*3e8*8.85418782e-12*1*cos(theta)*Emag_lambda.^2;
end

function [mag_xyz_lambda] = magnitudeE(monitor_name,variable)
    field = monitor_name.(variable);
    fieldX = field(:,1,:);
    fieldY = field(:,2,:);
    fieldZ = field(:,3,:);
    mag_xyz_lambda = sqrt(abs(fieldX).^2+abs(fieldY).^2+abs(fieldZ).^2);
end

function [E_spatial] = reshape_spatial(monitor_name,variable,field_component,wlVal)
%     component: x,y,z = 1,2,3
%     wlVal = index in the wavelength value array
    s1 = length(monitor_name.x);
    s2 = length(monitor_name.y);
    s3 = length(monitor_name.z);
    
    E = monitor_name.(variable)(:,field_component,wlVal);
    
    E_spatial = reshape(E,[s1 s2 s3]);
end