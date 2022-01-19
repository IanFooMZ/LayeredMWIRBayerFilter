clear; clc; close all;
%idcs = strfind(pwd,'\'); mydir = pwd; newdir = mydir(1:idcs(end)-1);
newdir = pwd;
load([newdir,'\','ang_range_superimposed_gauss_phi.mat']);

quadrant = 2;
legendVals = ["Blue","Green (x-pol.)","Red","Green (y-pol.)"];
saveFilename = ['th10_phiVals_',char(legendVals(quadrant+1)),''];

fig = figure; hold on;
for xi = [2.5:2.5:17.5]
    xline(xi,'-.','HandleVisibility','off','Color','#505050' ...
    ,'LineWidth', 1.0);
end
realBool = 0;

fig = addSpec(fig, th0phi0, quadrant, 0, '');
% fig = addSpec(fig, th5phi0, quadrant, -1, '');
% fig = addSpec(fig, th10phi0, quadrant, -2, '');
% fig = addSpec(fig, th15phi0, quadrant, 1, '');

% fig = addSpec(fig, th5phi0, quadrant, -1, '30x30um, f30');
% % fig = addSpec(fig, th5l40f40, quadrant, 0, '40x40um, f40');
% % fig = addSpec(fig, th5l40f30, quadrant, 0, '40x40um, f30');
% fig = addSpec(fig, th5l50f50, quadrant, 0, '50x50um, f50');
% fig = addSpec(fig, th5l50f20, quadrant, 1, '50x50um, f20');
% fig = addSpec(fig, th5l50f30, quadrant, 0, '50x50um, f30');

% fig = addSpec(fig,th0data,quadrant,0,'');
% fig = addSpec(fig,th0phi0,quadrant,0,'');
% fig = addSpec(fig,th5phi0,quadrant,0,'; \phi = 0');
% fig = addSpec(fig,th5phi30,quadrant,0,'; \phi = 30');
% fig = addSpec(fig,th5phi45,quadrant,0,'; \phi = 45');
% fig = addSpec(fig,th5phi60,quadrant,0,'; \phi = 60');
% fig = addSpec(fig,th5phi90,quadrant,0,'; \phi = 90');
% fig = addSpec(fig,th5phi180,quadrant,0,'; \phi = 180');
% fig = addSpec(fig,th5phi225,quadrant,0,'; \phi = 225');
% fig = addSpec(fig,th5phim30,quadrant,0,'; \phi = -30');
% fig = addSpec(fig,th5phim45,quadrant,0,'; \phi = -45');
% fig = addSpec(fig,th5phim90,quadrant,0,'; \phi = -90');
fig = addSpec(fig,th10phi0,quadrant,-2,'; \phi = 0');
fig = addSpec(fig,th10phi45,quadrant,1,'; \phi = 45');
fig = addSpec(fig,th10phi90,quadrant,2,'; \phi = 90');
fig = addSpec(fig,th10phi135,quadrant,0,'; \phi = 135');
fig = addSpec(fig,th10phi180,quadrant,0,'; \phi = 180');
fig = addSpec(fig,th10phi225,quadrant,0,'; \phi = 225');
fig = addSpec(fig,th10phim45,quadrant,0,'; \phi = -45');
fig = addSpec(fig,th10phim90,quadrant,2,'; \phi = -90');

xlim([0 17.5]); xticks([0:5:20]); %ylim([0.15 0.45]);
xlabel('\theta'); ylabel('Sorting Efficiency');
legend = legend('Location', 'southwest');
title(['Angular Range: ',char(legendVals(quadrant+1))]);
% set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
set(gcf,'position',[0 0 1920 1440]);

lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
  lines(i).MarkerSize = 2.0;
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)

realStr = '';
if realBool == true
    realStr = '_real';
end

exportgraphics(gca,[saveFilename,realStr,'.png']);
savefig([saveFilename,realStr,'.fig']);



function fig = addSpec(fig,struct,quadrant,disp,addOnStr)
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
    
    if contains(addOnStr,'40')
        specVal = 1.0*specVal;
    end
%     if struct.thetaOrig == 15
%         specVal = 1.08*specVal;
%     end
    
    plt = plot(struct.thetaVals+disp*h, specVal ...
        , 'o-', 'DisplayName', ['Optimized for ',num2str(struct.thetaOrig),'°',addOnStr] ...
        , 'MarkerSize',2);
%     if struct.thetaOrig == 0
%         plt.Color = '#BDB4AE';
%         plt.LineStyle = '-.';
%     end
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