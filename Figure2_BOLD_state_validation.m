clc;clear
close all
codepath = 'H:\ECoG\code_Trange';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\';

Excel = fullfile(WholePath,'CC_State_Validation_hipp_cortex.xlsx');
[~,~,CellData1] = xlsread(Excel,'Cortex_NREM');
[~,~,CellData2] = xlsread(Excel,'Cortex_REM');
F = figure('Position', [680 549 522 446]);
for lp=1:5   
data = cell2mat(CellData1(2:end,lp));
data(isnan(data))=[];
[h, ~] = raincloud_plot(data, 0.05+lp/3,'color',[17,140,56]/255, ...
    'box_on',1,'alpha',1);
delete(h{1});
hold on;
ylim([-0.35 0.35])
data = cell2mat(CellData2(2:end,lp));
data(isnan(data))=[];
[h, ~] = raincloud_plot(data, 0.15+lp/3,'color',[9,102,164]/255, ...
    'box_on',1,'alpha',1);
delete(h{1});

end
ylim([0 1.7]);xlim([-6 8])
view(90,-90)
set(gca,'tickdir','out','ticklength',[0.04 1],'linewidth',2,'fontsize',19);
saveas(F,fullfile(WholePath,'Activation','Validataion_Cortex_ECoG.emf'));

Excel = fullfile(WholePath,'CC_State_Validation_hipp_cortex.xlsx');
[~,~,CellData1] = xlsread(Excel,'HIPP_NREM');
[~,~,CellData2] = xlsread(Excel,'HIPP_REM');
F = figure('Position', [680 549 522 446]);
for lp=1:5   
data = cell2mat(CellData1(2:end,lp));
data(isnan(data))=[];
[h, ~] = raincloud_plot(data, 0.05+lp/3,'color',[17,140,56]/255, ...
    'box_on',1,'alpha',1);
delete(h{1});
hold on;
ylim([-0.35 0.35])

data = cell2mat(CellData2(2:end,lp));
data(isnan(data))=[];
[h, ~] = raincloud_plot(data, 0.15+lp/3,'color',[9,102,164]/255, ...
    'box_on',1,'alpha',1);
delete(h{1});

end
ylim([0 1.7]);xlim([-6 8])
view(90,-90)
set(gca,'tickdir','out','ticklength',[0.04 1],'linewidth',2,'fontsize',19);
saveas(F,fullfile(WholePath,'Activation','Validataion_HIPP_LFP.emf'));


Excel = fullfile(WholePath,'CC_State_Validation_hipp_cortex.xlsx');
[~,~,CellData] = xlsread(Excel,'State_CC');
F = figure;
for lp=1:5
plot([0 CellData{lp+2,2}],[1 1]*(-lp),'k-','linewidth',2);
hold on;
if CellData{lp+2,3}>0.05
    plot([CellData{lp+2,2}],[1]*(-lp),'o','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[150,150,150]/255);
else    
    plot([CellData{lp+2,2}],[1]*(-lp),'o','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[17,140,56]/255);
end

plot([0 CellData{lp+2,4}],[1 1]*(-lp-0.3),'k-','linewidth',2);
hold on;
if CellData{lp+2,5}>0.05
    plot([CellData{lp+2,4}],[1]*(-lp-0.3),'gs','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[150,150,150]/255);
else
    plot([CellData{lp+2,4}],[1]*(-lp-0.3),'gs','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[9,102,164]/255);
end
end
plot([0 0],[-6 0],'k','linewidth',2)
ylim([-6 0])
xlim([-0.5 0.5])
set(gca,'ytick','','box','off','tickdir','out',...
    'ticklength',[0.04 1],'linewidth',2,'fontsize',19)
saveas(F,fullfile(WholePath,'Activation','Validataion_HIPP_State_CC.emf'));


Excel = fullfile(WholePath,'CC_State_Validation_hipp_cortex.xlsx');
[~,~,CellData] = xlsread(Excel,'State_CC');
F = figure;
for lp=1:5
plot([0 CellData{lp+2,6}],[1 1]*(-lp),'k-','linewidth',2);
hold on;
if CellData{lp+2,7}>0.05
    plot([CellData{lp+2,6}],[1]*(-lp),'o','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[150,150,150]/255);
else    
    plot([CellData{lp+2,6}],[1]*(-lp),'o','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[17,140,56]/255);
end

plot([0 CellData{lp+2,8}],[1 1]*(-lp-0.3),'k-','linewidth',2);
hold on;
if CellData{lp+2,9}>0.05
    plot([CellData{lp+2,8}],[1]*(-lp-0.3),'gs','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[150,150,150]/255);
else
    plot([CellData{lp+2,8}],[1]*(-lp-0.3),'gs','MarkerSize',15,...
        'MarkerEdgeColor','none','MarkerFaceColor',[9,102,164]/255);
end
end
plot([0 0],[-6 0],'k','linewidth',2)
ylim([-6 0])
xlim([-0.5 0.5])
set(gca,'ytick','','box','off','tickdir','out',...
    'ticklength',[0.04 1],'linewidth',2,'fontsize',19)
saveas(F,fullfile(WholePath,'Activation','Validataion_Cortex_State_CC.emf'));

% state CC scatter


Excel = fullfile(WholePath,'CC_State_Validation_hipp_cortex.xlsx');
[~,~,CellData1] = xlsread(Excel,'Cortex_NREM');
[~,~,CellData2] = xlsread(Excel,'Cortex_REM');
F = figure('Position', [680 42 385 953]);
for lp=1:5   
    subplot(5,2,lp*2-1)
    data = cell2mat(CellData1(2:end,lp));
    data(isnan(data))=[];
    y = cell2mat(CellData1(2:end,6));
    y(isnan(y))=[];
    p = polyfit(data,y*100,1);
    y_ = polyval(p,[min(data) max(data)]);
    plot(data,y*100,'o','MarkerEdgeColor','none','MarkerSize',3,'MarkerFaceColor',[17,140,56]/255);
    hold on;plot([min(data) max(data)],y_,'k-','linewidth',1.5);
    ylim([-1 2]*1.5);%xlim([-3 9])
    [CC,p0] = corr(data,y);
    yli=ylim;xli=xlim;
    text(xli(1),yli(2),{['C.C.=',num2str(CC,'%.02f')];['p=',num2str(p0,'%.04f')]},'color','red');
    set(gca,'tickdir','in','ticklength',[0.03 1],'linewidth',1,'fontsize',15);

    subplot(5,2,lp*2-0)
    data = cell2mat(CellData2(2:end,lp));
    data(isnan(data))=[];
    y = cell2mat(CellData2(2:end,6));
    y(isnan(y))=[];
    p = polyfit(data,y*100,1);
    y_ = polyval(p,[min(data) max(data)]);
    plot(data,y*100,'o','MarkerEdgeColor','none','MarkerSize',3,'MarkerFaceColor',[9,102,164]/255);
    hold on;plot([min(data) max(data)],y_,'k-','linewidth',1.5);
    ylim([-1 2]*1.5);%xlim([-3 9])
    [CC,p0] = corr(data,y);
    yli=ylim;xli=xlim;
    text(xli(1),yli(2),{['C.C.=',num2str(CC,'%.02f')];['p=',num2str(p0,'%.04f')]},'color','red');
    set(gca,'tickdir','in','ticklength',[0.03 1],'linewidth',1,'fontsize',15);
end
saveas(F,fullfile(WholePath,'Activation','Validataion_Cortex_BOLD_Scatter.emf'));




Excel = fullfile(WholePath,'CC_State_Validation_hipp_cortex.xlsx');
[~,~,CellData1] = xlsread(Excel,'HIPP_NREM');
[~,~,CellData2] = xlsread(Excel,'HIPP_REM');
F = figure('Position', [680 42 385 953]);
for lp=1:5   
    subplot(5,2,lp*2-1)
    data = cell2mat(CellData1(2:end,lp));
    data(isnan(data))=[];
    y = cell2mat(CellData1(2:end,6));
    y(isnan(y))=[];
    p = polyfit(data,y*100,1);
    y_ = polyval(p,[min(data) max(data)]);
    plot(data,y*100,'o','MarkerEdgeColor','none','MarkerSize',3,'MarkerFaceColor',[17,140,56]/255);
    hold on;plot([min(data) max(data)],y_,'k-','linewidth',1.5);
    ylim([-1 2]*1.5);%xlim([-3 9])
    [CC,p0] = corr(data,y);
    yli=ylim;xli=xlim;
    text(xli(1),yli(2),{['C.C.=',num2str(CC,'%.02f')];['p=',num2str(p0,'%.04f')]},'color','red');
    set(gca,'tickdir','in','ticklength',[0.03 1],'linewidth',1,'fontsize',15);

    subplot(5,2,lp*2-0)
    data = cell2mat(CellData2(2:end,lp));
    data(isnan(data))=[];
    y = cell2mat(CellData2(2:end,6));
    y(isnan(y))=[];
    p = polyfit(data,y*100,1);
    y_ = polyval(p,[min(data) max(data)]);
    plot(data,y*100,'o','MarkerEdgeColor','none','MarkerSize',3,'MarkerFaceColor',[9,102,164]/255);
    hold on;plot([min(data) max(data)],y_,'k-','linewidth',1.5);
    ylim([-1 2]*1.5);%xlim([-3 9])
    [CC,p0] = corr(data,y);
    yli=ylim;xli=xlim;
    text(xli(1),yli(2),{['C.C.=',num2str(CC,'%.02f')];['p=',num2str(p0,'%.04f')]},'color','red');
    set(gca,'tickdir','in','ticklength',[0.03 1],'linewidth',1,'fontsize',15);
end
saveas(F,fullfile(WholePath,'Activation','Validataion_HIPP_BOLD_Scatter.emf'));
