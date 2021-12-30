%% Calculates sorting efficiency for each quadrant
% Definition: Fraction of incident power reaching target quadrant
% References: 
% https://support.lumerical.com/hc/en-us/articles/360034409554
% https://kx.lumerical.com/t/transforming-datasets-as-structures-to-matlab/2576/5
clear; clc; close all;
thetaVals = [36:41];%[0:2:16];

for k = [1:length(thetaVals)]
    theta = thetaVals(k);

    filename = 'sort_spec_bfast_data_z';
    load([filename, num2str(theta),'.mat']);

    wlVals = 1e6*E_fm0.lambda;

    %% Overall Incident Power
    Emag_fp0 = integrateOverSpace(E_fp0,'E',theta);

    %% Incident Power for Each Quadrant
    Emag_tm0 = integrateOverSpace(E_tm0,'E',theta);
    Emag_tm1 = integrateOverSpace(E_tm1,'E',theta);
    Emag_tm2 = integrateOverSpace(E_tm2,'E',theta);
    Emag_tm3 = integrateOverSpace(E_tm3,'E',theta);

    %% E-Field at Each Focal Monitor
    Emag_fm0 = integrateOverSpace(E_fm0,'E',theta);
    Emag_fm1 = integrateOverSpace(E_fm1,'E',theta);
    Emag_fm2 = integrateOverSpace(E_fm2,'E',theta);
    Emag_fm3 = integrateOverSpace(E_fm3,'E',theta);

    %% Plot Sorting Spectrum
    fig = figure; hold on;
    plot(wlVals,Emag_tm0./Emag_fp0,'b','DisplayName','Blue');
    plot(wlVals,Emag_tm1./Emag_fp0,'g','DisplayName','Green,x-pol');
    plot(wlVals,Emag_tm2./Emag_fp0,'r','DisplayName','Red');
    plot(wlVals,Emag_tm3./Emag_fp0,'Color',1/255*[40,94,25],'DisplayName','Green,y-pol');

    xlabel('Wavelength (um)');
    ylabel('Sorting Efficiency');
%     ylim([0.05,0.55]);
    ylim([0.1,0.55]);
    leg = legend('Location', 'north');
    title(['Spectrum in Each Quadrant, ',num2str(theta),'° incidence']);
    set(gcf,'position',[361.0000  226.3333  675.3333  392.6667]);
    saveas(fig,['sort_spec_bfast_theta0_z', num2str(theta), '.png']);
    close all;

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