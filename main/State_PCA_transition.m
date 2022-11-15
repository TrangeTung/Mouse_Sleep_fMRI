close all
clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\pca\';


load(fullfile('H:\ECoG\dcc_yyl\','states_1s.mat'));


StateName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

% phase 
GS_S1=[]; GS_S2=[]; 
GS_S3=[]; GS_S4=[];
St_S1=[]; St_S2=[]; 
St_S3=[]; St_S4=[];
CutBins = 30;
CutBins_2nd = 10;

for idx=1:46
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));

    %xState = ALL(:,idx);
    xState =  all_state;
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    Ss = State;
    
    dS = Ss-[Ss(1);Ss(1:end-1)];
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(dS)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.21);
    beta = filtfilt(b,a,beta')';
    
%     y = hilbert(beta);
%     phase = angle(y);
%     beta = phase;
    
    T = find(dS==+2); % state 1
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==1&unique(Ss(T(tl)+(1:CutBins_2nd)))==3
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S1 = cat(3,GS_S1,Gc(:,:,mean(St,1)~=0));
    St_S1 = cat(2,St_S1,St(:,mean(St,1)~=0));

    T = find(dS==-2); % state 2
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S2 = cat(3,GS_S2,Gc(:,:,mean(St,1)~=0));
    St_S2 = cat(2,St_S2,St(:,mean(St,1)~=0));
    
    
    T = find(dS==+1); % state 3
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==4
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S3 = cat(3,GS_S3,Gc(:,:,mean(St,1)~=0));
    St_S3 = cat(2,St_S3,St(:,mean(St,1)~=0));

    
    T = find(dS==-3); % state 4
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==4&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S4 = cat(3,GS_S4,Gc(:,:,mean(St,1)~=0));
    St_S4 = cat(2,St_S4,St(:,mean(St,1)~=0));

    
end

StateName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};
PC_ColorVec = {[217,83,25];[0 ,104,60];[140,182,44];[171,131,71]};



for ss=1:4
    
    x_ = (-CutBins:CutBins)*2;

    if ss==1;Bs=[1 3];end
    if ss==2;Bs=[3 1];end
    if ss==3;Bs=[3 4];end
    if ss==4;Bs=[4 1];end
    
    eval(['GS=GS_S',num2str(ss),';']);
    eval(['St=St_S',num2str(ss),';']);
    
    
    F = figure( 'Position', [680 466 469 512]);
    subplot(5,1,1);
    y_ = [mean(St==Bs(1),2),mean(St==Bs(2),2)];
    bar(x_,y_,'stacked','edgecolor','none');hold on;
    xlim([min(x_) max(x_)]);
    
    subplot(5,1,2:5);
    for pl=1:4
        S = squeeze(GS(pl,:,:));
        
        Ms = mean(S,2);
        SEM = std(S,0,2)./sqrt(size(S,2));
        
        plot(x_,Ms,'color',PC_ColorVec{pl}/255,'linewidth',2);
        hold on;
        patch([x_(:);flipud(x_(:))],[Ms+SEM;flipud(Ms-SEM)],...
            ones([numel(x_)*2,1]),...
            'facecolor',PC_ColorVec{pl}/255,...
            'edgecolor','none','facealpha',.25);
        
    end
    xlim([min(x_) max(x_)]);
    set(gca,'tickdir','out');
    saveas(F,fullfile(WholePath,['State_Trans_PCweight_',StateName{ss},'.emf']));
end



PC_ColorVec = {[224,102,47];[19,103,54];[134,183,35];[172,125,67]};

F = figure( 'Position', [280 066 765 641]);

for ss=1:4
    
    x_ = (-CutBins:CutBins);

    if ss==1;Bs=[1 3];end
    if ss==2;Bs=[3 1];end
    if ss==3;Bs=[3 4];end
    if ss==4;Bs=[4 1];end
    
    eval(['GS=GS_S',num2str(ss),';']);
    eval(['St=St_S',num2str(ss),';']);
    
    y_ = - mean(St==Bs(1),2) + mean(St==Bs(2),2);
    y_(31) = y_(1)/2+y_(end)/2;

    for pl=1:4
        
        
        S = squeeze(GS(pl,:,:));
        
        subplot(4,4,(pl-1)*4+ss);
        polarhistogram(S(31,:),20);
    end
end




% shift CC


PC_ColorVec = {[224,102,47];[19,103,54];[134,183,35];[172,125,67]};

F = figure( 'Position', [680 466 765 641]);

for ss=1:4
    
    x_ = (-CutBins:CutBins);

    if ss==1;Bs=[1 3];end
    if ss==2;Bs=[3 1];end
    if ss==3;Bs=[3 4];end
    if ss==4;Bs=[4 1];end
    
    eval(['GS=GS_S',num2str(ss),';']);
    eval(['St=St_S',num2str(ss),';']);
    
    y_ = - mean(St==Bs(1),2) + mean(St==Bs(2),2);
    y_(31) = y_(1)/2+y_(end)/2;
  
    

    for pl=1:4
        
        
        S = squeeze(GS(pl,:,:));
        
        Ms = mean(S,2);
        K(:,(ss-1)*4+pl)=Ms;
        
        CC = zeros(numel(x_),1);
        for sl=1:numel(x_)
           sx = circshift(y_,x_(sl),1);
           CC(sl) = corr(sx(:),Ms(:));
        end
        
        subplot(4,4,(pl-1)*4+ss);
        plot(x_*2,CC,'color',PC_ColorVec{pl}/255,'linewidth',2);
        
        hold on;plot([-100 100],[0 0]);
        hold on;plot([0 0],[-1 1]);
        xlim([min(x_) max(x_)]*2.1);
        ylim([-1 1])
        set(gca,'tickdir','out');

        hold on;
        %CC([1:10,end-10:end])=0;
        MNa = x_(find(CC==max(CC)))*2;
        MNi = x_(find(CC==min(CC)))*2;
        bb=bar([MNa],[max(CC)],5,'r','edgecolor','none');
        xtips1 = bb(1).XEndPoints;
        ytips1 = bb(1).YEndPoints;
        labels1 = string(bb(1).XData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        bb=bar([MNi],[min(CC)],5,'b','edgecolor','none');
        xtips1 = bb(1).XEndPoints;
        ytips1 = bb(1).YEndPoints;
        labels1 = string(bb(1).XData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','top');

        title(['State ', num2str(ss),' PC ',num2str(pl)]);
        if ss~=1|pl~=1; box off;axis off;end
    end
end

saveas(F,fullfile(WholePath,['State_Trans_PCweight.emf']));
%}