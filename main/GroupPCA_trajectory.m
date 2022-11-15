clc;clear
codepath = 'F:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'F:\ECoG\pca\';

close all

for loop=1:4

    if loop==1;XYZ=[1 2 3];end
    if loop==2;XYZ=[1 2 4];end
    if loop==3;XYZ=[1 3 4];end
    if loop==4;XYZ=[2 3 4];end

%% averaged
F = figure('Position', [680 478 554 500]);

State = {'awake';'nrem';'rem'};
ColorVec = {[215,218,149];[142,212,156];[139,148,195]};
for ss=[1 2 3]
    
    dest = fullfile(WholePath,'Trajectory');
    load(fullfile(dest,['Trajectory_',State{ss},'.mat']));
    
    S = nanmean(St,3);
    scatter3(S(XYZ(1),:)',S(XYZ(2),:)',S(XYZ(3),:)',25,'MarkerFacealpha',0.9,...
        'MarkerEdgeColor','none','MarkerFaceColor',ColorVec{ss}/255);
    hold on;
    %plot3(S(XYZ(1),:)',S(XYZ(2),:)',S(XYZ(3),:)','k-');
    
end

State = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};
ColorMap{1} = [linspace(215,142,202)',linspace(218,212,202)',linspace(149,156,202)'];
ColorMap{2} = [linspace(142,215,202)',linspace(212,218,202)',linspace(156,149,202)'];
ColorMap{3} = [linspace(142,139,202)',linspace(212,148,202)',linspace(156,195,202)'];
ColorMap{4} = [linspace(139,215,202)',linspace(148,218,202)',linspace(195,149,202)'];

for ss=[1 2 3 4]
    
    dest = fullfile(WholePath,'Trajectory');
    load(fullfile(dest,['Trajectory_',State{ss},'.mat']));
    S = nanmean(St,3);
    scatter3(S(XYZ(1),[1:end,1])',S(XYZ(2),[1:end,1])',S(XYZ(3),[1:end,1])',25,ColorMap{ss}/255,...
        'filled','MarkerFacealpha',0.9,...
        'MarkerEdgeColor','none');
    hold on;
    %plot3(S(XYZ(1),:)',S(XYZ(2),:)',S(XYZ(3),:)','k-');
end

xlabel(['PC',num2str(XYZ(1))]);
ylabel(['PC',num2str(XYZ(2))]);
zlabel(['PC',num2str(XYZ(3))]);
set(gca,'linewidth',1,'fontsize',15);
view([200 20]);

grid on;

    if loop==1;view([040 20]);end
    if loop==2;view([200 20]);end
    if loop==3;view([223 47]);end
    if loop==4;view([300 50]);end




%% scatter
F = figure('Position', [680 478 554 500]);

State = {'awake';'nrem';'rem'};
ColorVec = {[215,218,149];[142,212,156];[139,148,195]};
for ss=[1 2 3]
    
    dest = fullfile(WholePath,'Trajectory');
    load(fullfile(dest,['Trajectory_',State{ss},'.mat']));
    bs = 10;
    if ss==3;bs=10;end
    NUM = floor(size(St,3)/bs);
    for bl=1:NUM
        Cs = (bl-1)*bs+1:bl*bs;
        S = nanmean(St(:,:,Cs),3);
        scatter3(S(XYZ(1),:)',S(XYZ(2),:)',S(XYZ(3),:)',25,'MarkerFacealpha',0.1,...
            'MarkerEdgeColor','none','MarkerFaceColor',ColorVec{ss}/255);
        hold on;
    end
    
end

State = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};
ColorMap{1} = [linspace(215,142,202)',linspace(218,212,202)',linspace(149,156,202)'];
ColorMap{2} = [linspace(142,215,202)',linspace(212,218,202)',linspace(156,149,202)'];
ColorMap{3} = [linspace(142,139,202)',linspace(212,148,202)',linspace(156,195,202)'];
ColorMap{4} = [linspace(139,215,202)',linspace(148,218,202)',linspace(195,149,202)'];

for ss=[1 2 3 4]
    
    dest = fullfile(WholePath,'Trajectory');
    load(fullfile(dest,['Trajectory_',State{ss},'.mat']));
    if ss>=3; bs = 10;end
    if ss<=2; bs = 10;end
    NUM = floor(size(St,3)/bs);
    for bl=1:NUM
        Cs = (bl-1)*bs+1:bl*bs;
        S = nanmean(St(:,:,Cs),3);
        scatter3(S(XYZ(1),[1:end,1])',S(XYZ(2),[1:end,1])',S(XYZ(3),[1:end,1])',25,ColorMap{ss}/255,...
            'filled','MarkerFacealpha',0.1,...
            'MarkerEdgeColor','none');
        hold on;
    end
    
end

xlabel(['PC',num2str(XYZ(1))]);
ylabel(['PC',num2str(XYZ(2))]);
zlabel(['PC',num2str(XYZ(3))]);
set(gca,'linewidth',1,'fontsize',15);
view([40 20]);


    if loop==1;view([040 20]);end
    if loop==2;view([200 20]);end
    if loop==3;view([223 47]);end
    if loop==4;view([300 50]);end




% saveas(F,fullfile(dest,['Trajectory_full_',num2str(XYZ),'.tiff']));


end