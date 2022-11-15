close all
clc;clear
codepath = 'F:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'F:\ECoG\pca\';


%load(fullfile('F:\ECoG\dcc_yyl\','states_tr.mat'));
load(fullfile('F:\ECoG\dcc_yyl\','states_1s.mat'));

GS_S1=[]; GS_S2=[]; 
GS_S3=[];


for idx=1:46
    xState = ALL(:,idx);
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    
    dS = State-[State(1);State(1:end-1)];
    
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(State)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.2);
    %beta = filtfilt(b,a,beta')';
    
    X=(beta(:,State==1)); NUM = floor(size(X,2)/100); Y=[];
    for nl=1:NUM;Y(:,nl)=mean(X(:,(nl-1)*100+(1:100)),2);end
    GS_S1 = cat(2,GS_S1,Y);
    
    X=(beta(:,State==3)); NUM = floor(size(X,2)/100); Y=[];
    for nl=1:NUM;Y(:,nl)=mean(X(:,(nl-1)*100+(1:100)),2);end
    GS_S2 = cat(2,GS_S2,Y);

    X=(beta(:,State==4)); NUM = floor(size(X,2)/100); Y=[];
    for nl=1:NUM;Y(:,nl)=mean(X(:,(nl-1)*100+(1:100)),2);end
    GS_S3 = cat(2,GS_S3,Y);


end

clear data

F = figure('Position', [1275 573 388 352]);
for loop=1:3
    subplot(1,3,loop);
    eval(['X=GS_S',num2str(loop),';']);
    h = raincloud_plot(X(1,:), 'box_on', 1, 'color', [0.5 0.5 0.5]);
    xlim([-3 3]);view(-90,90);
    hold on;
end
% saveas(F,fullfile(WholePath,'State_PC1.emf'));

F = figure('Position', [1275 573 388 352]);
for loop=1:3
    subplot(1,3,loop);
    eval(['X=GS_S',num2str(loop),';']);
    h = raincloud_plot(X(2,:), 'box_on', 1, 'color', [0.5 0.5 0.5]);
    xlim([-3 3]);view(-90,90);
    hold on;
end
% saveas(F,fullfile(WholePath,'State_PC2.emf'));

F = figure('Position', [1275 573 388 352]);
for loop=1:3
    subplot(1,3,loop);
    eval(['X=GS_S',num2str(loop),';']);
    h = raincloud_plot(X(3,:), 'box_on', 1, 'color', [0.5 0.5 0.5]);
    xlim([-3 3]);view(-90,90);
    hold on;
end
% saveas(F,fullfile(WholePath,'State_PC3.emf'));

F = figure('Position', [1275 573 388 352]);
for loop=1:3
    subplot(1,3,loop);
    eval(['X=GS_S',num2str(loop),';']);
    h = raincloud_plot(X(4,:), 'box_on', 1, 'color', [0.5 0.5 0.5]);
    xlim([-3 3]);view(-90,90);
    hold on;
end
% saveas(F,fullfile(WholePath,'State_PC4.emf'));


