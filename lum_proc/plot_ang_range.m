%% Calculates sorting efficiency for each quadrant
% Definition: Fraction of incident power reaching target quadrant
% References: 
% https://support.lumerical.com/hc/en-us/articles/360034409554
% https://kx.lumerical.com/t/transforming-datasets-as-structures-to-matlab/2576/5
clear; clc; close all;
thetaOrig = 5;
% thetaVals = [thetaOrig-0:1.25:thetaOrig+20];
thetaVals = [thetaOrig-15:1.25:thetaOrig+15];
peakInd = find(thetaVals==thetaOrig);
sourceType = 'gauss';
dp = 0;

file{1} = [pwd,'\'];%,num2str(thetaOrig,'%.1f'),'_inverse_design\'];
file{2} = 'sortspecdata';
file{4} = ['optang',num2str(thetaOrig,'%.1f'),'_th'];
file{3} = [sourceType];%,'inSi'];%, '_xpol'];
fn = formFileName(file,thetaOrig,dp);

try
    load([fn,'.mat']);
catch ME
    load([fn(1:end-length(file{4})-1),'.mat']);
end
wlVals = 1e6*E_fm0.lambda;

Emag_fp0 = zeros(length(thetaVals), length(wlVals));
Emag_tm0 = Emag_fp0;
Emag_tm1 = Emag_fp0;
Emag_tm2 = Emag_fp0;
Emag_tm3 = Emag_fp0;

for k = 1:length(thetaVals)
    theta = thetaVals(k);
    if theta<0
        fn = formFileName(file,ceil(theta),dp);
    else
        fn = formFileName(file,floor(theta),dp);
    end
    try
        load([fn,'.mat']);
    catch ME
        load([fn(1:end-length(file{4})-1),'.mat']);
    end
    wlVals = 1e6*E_fm0.lambda;

    %% Overall Incident Power
    Emag_fp0(k,:) = integrateOverSpace(E_fp0,'E',theta);

    %% Incident Power for Each Quadrant
    Emag_tm0(k,:) = integrateOverSpace(E_tm0,'E',theta);
    Emag_tm1(k,:) = integrateOverSpace(E_tm1,'E',theta);
    Emag_tm2(k,:) = integrateOverSpace(E_tm2,'E',theta);
    Emag_tm3(k,:) = integrateOverSpace(E_tm3,'E',theta);

%     %% E-Field at Each Focal Monitor
%     Emag_fm0 = integrateOverSpace(E_fm0,'E',theta);
%     Emag_fm1 = integrateOverSpace(E_fm1,'E',theta);
%     Emag_fm2 = integrateOverSpace(E_fm2,'E',theta);
%     Emag_fm3 = integrateOverSpace(E_fm3,'E',theta);

end


%% Normalize by focal plane, then take Maximum along wavelength dimension
%% (designed peak), then normalize 
% Emag_tm0 = max(Emag_tm0./Emag_fp0,[],2);
% Emag_tm1 = max(Emag_tm1./Emag_fp0,[],2);
% Emag_tm2 = max(Emag_tm2./Emag_fp0,[],2);
% Emag_tm3 = max(Emag_tm3./Emag_fp0,[],2);
Emag_tm0 = findMax(Emag_tm0./Emag_fp0, peakInd);
Emag_tm1 = findMax(Emag_tm1./Emag_fp0, peakInd);
Emag_tm2 = findMax(Emag_tm2./Emag_fp0, peakInd);
Emag_tm3 = findMax(Emag_tm3./Emag_fp0, peakInd);

% Save to struct
saveData = struct;
saveData.wlVals = wlVals;
saveData.peakInd = peakInd;
saveData.thetaVals = thetaVals;
saveData.thetaOrig = thetaOrig;
saveData.Emag_fp0 = Emag_fp0;
saveData.Emag_tm0 = Emag_tm0;
saveData.Emag_tm1 = Emag_tm1;
saveData.Emag_tm2 = Emag_tm2;
saveData.Emag_tm3 = Emag_tm3;

idcs = strfind(pwd,'\'); mydir = pwd; newdir = mydir(1:idcs(end)-1);
th5len40f30data = saveData; save([newdir,'\','ang_range_superimposed_',sourceType,'_Footprint.mat'] ...
    ,'th5len40f30data','-append');

Emag_tm0 = Emag_tm0./Emag_tm0(peakInd);
Emag_tm1 = Emag_tm1./Emag_tm1(peakInd);
Emag_tm2 = Emag_tm2./Emag_tm2(peakInd);
Emag_tm3 = Emag_tm3./Emag_tm3(peakInd);

%% Make it Look Good
dth = thetaVals(2)-thetaVals(1);

Emag_tm0 = Emag_tm0./max(Emag_tm0);
Emag_tm1 = Emag_tm1./max(Emag_tm1);
Emag_tm2 = Emag_tm2./max(Emag_tm2);

%% Plot Sorting Spectrum
fig = figure; hold on;
xline(thetaOrig,'-','HandleVisibility','off','Color','#505050' ...
    ,'LineWidth', 2.0);
yline(0.5,'--','HandleVisibility','off','Color','#505050' ...
    ,'LineWidth', 2.0);

intensity = 230;
plot(thetaVals-0*dth,Emag_tm0,'o-','Color',[0 0 intensity]./255,'DisplayName','Blue');
plot(thetaVals-1*dth,Emag_tm1,'o-','Color',[0 intensity 0]./255,'DisplayName','Green, x-pol');
plot(thetaVals+0*dth,Emag_tm2,'o-','Color',[intensity 0 0]./255,'DisplayName','Red');
% plot(thetaVals,Emag_tm3,'Color',1/255*[40,94,25],'DisplayName','Green,y-pol');

xlabel('Angle of Incidence (°)');
ylabel('Sorting Efficiency (Normalized)');
xlim([thetaVals(2) thetaVals(end-1)]); 
ylim([0.4,1.05]);
legend = legend('Location', 'northwest');
title(['Angular Range: Optimized at ',num2str(thetaOrig,"%.1f"),'°']);

lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
  lines(i).MarkerSize = 2.0;
end
set(findall(gcf,'-property','FontSize'),'FontSize',16)

%set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
set(gcf,'position',[0 0 1920 1440]);
exportgraphics(gca,['angrange_th',num2str(thetaOrig),'.png']);

%% Functions
function fn = formFileName(file,theta,dp)
    switch dp
        case 0
            thetaStr = num2str(theta,'%.0f');
        case 1
            thetaStr = num2str(theta,'%.1f');
        case 2
            thetaStr = num2str(theta,'%.2f');
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