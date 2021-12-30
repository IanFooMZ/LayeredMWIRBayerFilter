clear; clc; close all;
load('ang_range_superimposed_tfsf.mat');

quadrant = 1;
legendVals = ["Blue","Green (x-pol.)","Red","Green (y-pol.)"];

fig = figure; hold on;
for xi = [2.5:2.5:10]
    xline(xi,'--','HandleVisibility','off','Color','#505050' ...
    ,'LineWidth', 2.0);
end
fig = addSpec(fig,th0data,quadrant,0);
fig = addSpec(fig,th2_5data,quadrant,-1);
fig = addSpec(fig,th5data,quadrant,-1);
fig = addSpec(fig,th7_5data,quadrant,-2);
fig = addSpec(fig,th10data,quadrant,-2); 

xlim([0 15]); %ylim([0.15 0.45]);
xlabel('\theta'); ylabel('Sorting Efficiency');
legend = legend('Location', 'northeast');
title(['Angular Range: ',char(legendVals(quadrant+1))]);
% set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
set(gcf,'position',[0 0 1920 1440]);

lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
  lines(i).MarkerSize = 2.0;
end
set(findall(gcf,'-property','FontSize'),'FontSize',16)

exportgraphics(gca,['angrangesuper_',char(legendVals(quadrant+1)),'.png']);

function fig = addSpec(fig,struct,quadrant,disp)
    fig = ancestor(fig,'figure');
    
    switch quadrant
        case 0
            specVal = struct.Emag_tm0;
        case 1
            specVal = struct.Emag_tm1;
        case 2
            specVal = struct.Emag_tm2;
        case 3
            specVal = struct.Emag_tm3;
    end
    
    h = struct.thetaVals(2) - struct.thetaVals(1);
    plot(struct.thetaVals+disp*h, specVal ...
        , 'o-', 'DisplayName', ['Optimized for ',num2str(struct.thetaOrig),'°'] ...
        , 'MarkerSize',2);
end
%% Functions
function [Emag_processed] = findMax(Emag,peakInd)
    [a, ind] = max(Emag(peakInd,:));
    Emag_processed = Emag(:,ind);
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