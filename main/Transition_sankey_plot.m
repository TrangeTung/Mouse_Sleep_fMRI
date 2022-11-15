clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\Activity_transition\';

State = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};
Excel = 'H:\ECoG\Submit\v3\tran_30s_raw.xlsx';


for sl=1:numel(State)
    [~,~,CellData] = xlsread(Excel,State{sl});
    D = CellData(2:end,3:end);
    D(strcmp(D,'NaN'))={NaN};
    Data = cell2mat(D);
    
    t = -30:2:30;
    
    F = figure('Position', [680 266 560 712]);
    subplot(3,2,[2 4 6]);
    gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
    drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
    defaultMap = [[0.9 0.9 0.9];gdmap;flipud(drmap);];
    imagesc(t,1:size(Data,1),Data);
    colormap(defaultMap); colorbar
    if sl==1|sl==2;caxis([-0.4 0.4]);end
    if sl==3|sl==4;caxis([-1.5 1.5]);end
    set(gca,'tickdir','out','fontsize',15);
    
    subplot(3,2,[1 3 5]);
    List = {};
    List(:,1) = CellData(2:end,2);
    List(:,2) = repmat({1},numel(1:size(Data,1)),1);
    List(:,3) = CellData(2:end,1);
    colorList = jet(size(Data,1));
    sankeyHdl=sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',List,'Color',colorList);
    axis off
    saveas(F,fullfile(WholePath,[State{sl},'_series1.emf']));
    close all
end

%% Global signal


StateName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

% phase 
GS_S1=[]; GS_S2=[]; 
GS_S3=[]; GS_S4=[];
CutBins = 15;
CutBins_2nd = 5;

for idx=1:46
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));

    %xState = ALL(:,idx);
    xState =  all_state;
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    Ss = State;
    
    dS = Ss-[Ss(1);Ss(1:end-1)];
    beta=load(fullfile('H:\ECoG\G_S\G_S',num2str(idx,'%02d'),'G_S.txt'));%
    beta(:,numel(dS)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.51);
%     beta = filtfilt(b,a,beta')';
    
    
    T = find(dS==+2); % state 1
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(1,CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==1&unique(Ss(T(tl)+(1:CutBins_2nd)))==3
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S1 = cat(3,GS_S1,Gc(:,:,:));

    T = find(dS==-2); % state 2
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(1,CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S2 = cat(3,GS_S2,Gc(:,:,:));
    
    
    T = find(dS==+1); % state 3
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(1,CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==4
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S3 = cat(3,GS_S3,Gc(:,:,:));

    
    T = find(dS==-3); % state 4
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(1,CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==4&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S4 = cat(3,GS_S4,Gc(:,:,:));

    
end


F = figure;
subplot(2,2,1);
plot(-30:2:30,squeeze(mean(GS_S1,3))); 
ylim([-.1 .1]*4); xlim([-30 30]);
hold on;plot([0  0],[-.1 .1]*8);plot([-30 30],[0 0]);
subplot(2,2,2);
plot(-30:2:30,squeeze(mean(GS_S2,3)));
ylim([-.1 .1]*4); xlim([-30 30]);
hold on;plot([0  0],[-.1 .1]*8);plot([-30 30],[0 0]);
subplot(2,2,3);
plot(-30:2:30,squeeze(mean(GS_S3,3)));
ylim([-.1 .1]*20); xlim([-30 30]);
hold on;plot([0  0],[-.1 .1]*20);plot([-30 30],[0 0]);
subplot(2,2,4);
plot(-30:2:30,squeeze(mean(GS_S4,3)));
ylim([-.1 .1]*20); xlim([-30 30]);
hold on;plot([0  0],[-.1 .1]*20);plot([-30 30],[0 0]);




filepath = 'H:\ECoG\Activity_transition\';
Temp3D = fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii');
for sl=1:numel(State)
    barlim = [-30 -0.1 0.1 30];
Func3D = fullfile(filepath,['X_',State{sl},'_timelag.nii']);
 Img_RGB=Colormap(   'statfile',Func3D,...
                    'bgmfile',Temp3D,...
                    'slice',208:-14:81,...
                    'bar_value',barlim,...
                    'denoi_profile',fullfile(codepath,'Colormap_3Dviewer','lmask_Mouse_X20.nii'),...
                    'dest',fullfile(filepath),...
                    'mapname',[State{sl},'_slice'],...
                    'cluster',10);
end