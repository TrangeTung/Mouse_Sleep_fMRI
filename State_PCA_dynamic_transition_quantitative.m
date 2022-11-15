close all
clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\pca\';


load(fullfile('H:\ECoG\dcc_yyl\','states_1s.mat'));


StateName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};



for shf = 00%:1000
% phase 
GS_S1=[]; GS_S2=[]; 
GS_S3=[]; GS_S4=[];
St_S1=[]; St_S2=[]; 
St_S3=[]; St_S4=[];
CutBins = 30;
CutBins_2nd = 10;
 
    

for idx=1:46
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    State =  all_state;
   
    State(State==2)=1;
    uniqueS = unique(State);
    
    dx = State-[State(1);State(1:end-1)];
    P = double(dx~=0);
    [~,LOCS] = findpeaks(P,'MinPeakHeight',0.5);
    LOCS=[1;LOCS;numel(P)+1];
    Ss={};
    for xl=1:numel(LOCS)-1
        Ss{1,xl}=State(LOCS(xl):LOCS(xl+1)-1);
    end
    if shf~=0;Ss=Ss(randperm(numel(Ss)));end
    xState = cat(1,Ss{:});
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    Ss = State;
    
    dS = Ss-[Ss(1);Ss(1:end-1)];
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(dS)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.21);
    beta = filtfilt(b,a,beta')';
    
    
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

dest = fullfile(WholePath,'StatetTransition_Quantitative');
if ~exist(dest,'dir');mkdir(dest);end
AVG_S1 = nanmean(GS_S1,3);
AVG_S2 = nanmean(GS_S2,3);
AVG_S3 = nanmean(GS_S3,3);
AVG_S4 = nanmean(GS_S4,3);

save(fullfile(dest,['Shuffle_',num2str(shf,'%04d'),'.mat']),'AVG_S1','AVG_S2','AVG_S3','AVG_S4');
end




clear Shuffle ;
for shf = 00:1000
    filepath = fullfile(WholePath,'StatetTransition_Quantitative');
    x = load(fullfile(filepath,['Shuffle_',num2str(shf,'%04d'),'.mat']));
    
    if shf==0; Raw = x; end
    if shf~=0; Shuffle(shf) = x; end
end

StateName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

close all
for st = [1,2,3,4];
    
    
    F = figure('Position', [442 714 1381 264]*1.2);
    for pc=1:4
        eval(['x1=Raw.AVG_S',num2str(st),'(',num2str(pc),',:,:);']);
        eval(['y1=cat(1,Shuffle.AVG_S',num2str(st),');']);
        y2 = y1(pc:100:end,:);
        
        subplot(1,4,pc);
        x = squeeze(x1);
        Ms = x;%mean(x,2); 
        t0 = 2*(-30:30)';
        if pc==1;ColorVec=[217,83,25]/255;end
        if pc==2;ColorVec=[0 ,104,60]/255;end
        if pc==3;ColorVec=[140,182,44]/255;end
        if pc==4;ColorVec=[171,131,71]/255;end
        p=plot(t0,Ms,'color',[ColorVec],'linewidth',2);
        hold on;
        patch('XData',[t0;flipud(t0)]','YData',...
            [mean(y2,1)-2.0*std(y2,0,1),fliplr(mean(y2,1)+2.0*std(y2,0,1))],...
            'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');
        set(gca,'box','off','tickdir','out','ticklength',[0.04 0.1])
        if st==1;ylim([-0.3 0.4]);end
        if st==2;ylim([-0.8 0.4]);end
        if st==3;ylim([-1.5 1.5]);end
        if st==4;ylim([-1.5 1.5]);end
    end
    saveas(F,fullfile('H:\ECoG\pca',['State_',StateName{st},'_control.emf']));
end



