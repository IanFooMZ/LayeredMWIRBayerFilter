%% Superimposes the angular ranges of respective runs where a different 
%% variable can be swept.
%% Can also be used to conduct a 2D sweep, where the 1st sweep parameter
%% has already been compiled into the file that this code loads;
%% This code will then present those plots superimposed on top of each other
clear; clc; close all;

%% Load the compile file
%idcs = strfind(pwd,'\'); mydir = pwd; newdir = mydir(1:idcs(end)-1);
newdir = pwd;
load([newdir,'\','angRange_suprImp_gauss_bRad.mat']);

%% Choose Appropriate Frequency Interval
lgdVals = ["Blue","Green (x-pol.)","Red","Green (y-pol.)"];

for quadrant = [0:2]

    %% Set up Figure Plotting
    realBool = 0; realStr = '';

    % Figure Properties
    fig = figure; ax = gca; hold on;
    xlim([0 17.5]); xticks([0:5:20]); %ylim([0.15 0.45]);
    xlabel('\theta'); ylabel('Sorting Efficiency');
    lgd = legend('Location', 'southwest');
    title(['Angular Range: ',char(lgdVals(quadrant+1))]);

    for xi = [2.5:2.5:17.5]
        xline(xi,'-.','HandleVisibility','off','Color','#505050' ...
        ,'LineWidth', 1.0);
    end

    if realBool
        realStr = '_real';
    end
    sfn = {['ars_',...
        char(lgdVals(quadrant+1)),...
        realStr]};

    %% Plot each Set of Lines
    fig = addSpec(fig, th5bRad_test, quadrant, realBool, 0, 'th5, bRad (test)');

    %% Figure Post-Processing
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
    close all;

end


function fig = addSpec(fig,struct,quadrant,realBool,disp,lgdStr)
    fig = ancestor(fig,'figure');

    specVal = struct.(['Parr_tm',num2str(quadrant,'%d')]);
    dswp = struct.sweepVals(2) - struct.sweepVals(1);

    if ~realBool
       % usually the input displacement is 0
       disp = disp + struct.ov(quadrant+1);
    end

    % Pre-Process
%     if contains(lgdStr,'135')
%         specVal = 1.05*specVal;
%     end

    % Plot actual line
    plt = plot(struct.sweepVals+disp*dswp, specVal ...
        , 'o-', 'DisplayName', ['Optimized for ',num2str(struct.thetaOrig),'°',lgdStr] ...
        , 'MarkerSize',2);

    % Post-Process
%     if struct.thetaOrig == 0
%         plt.Color = '#BDB4AE';
%         plt.LineStyle = '-.';
%     end
end