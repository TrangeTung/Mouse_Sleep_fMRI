
load(fullfile('H:\ECoG\dcc_yyl\','states_1s.mat'));


StateName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

% phase 
GS_S1=[]; GS_S2=[]; 
GS_S3=[]; GS_S4=[];
St_S1=[]; St_S2=[]; 
St_S3=[]; St_S4=[];

BOLD_ALL_S1 = zeros(100,31,46);
BOLD_ALL_S2 = zeros(100,31,46);
ECoG_ALL_S1 = zeros(2,61,46);
ECoG_ALL_S2 = zeros(2,61,46);
for idx=[1:15 17:46]
    %% BOLD
    CutBins = 30/2;
    CutBins_2nd = 10;

    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));

    %xState = ALL(:,idx);
    xState =  all_state;
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    Ss = State;
    
    WholePath = 'H:\ECoG\pca\';
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
    BOLD_ALL_S1(:,:,idx)=mean(Gc,3);
    
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
    BOLD_ALL_S2(:,:,idx)=mean(Gc,3);
    
    
    
    %% ECoG
    CutBins = 30;
    CutBins_2nd = 10*2;
    WholePath = 'H:\ECoG\yyl_power\';

    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    all_state(all_state==2)=1;
    cd(fullfile(WholePath,[num2str(idx)]));
    
    try
       load(fullfile(WholePath,[num2str(idx)],['x_S_100_5.mat']));
       S=S;t=t;f=f;
    catch
        load(fullfile(WholePath,[num2str(idx)],['x_S_100_7.mat'])); 
        S=S;t=t;f=f;
        1;
    end
    
    X=10*log10(S); X = abs(X);
    fc = 11*(1:8);
    for q=1:8;X(:,f<fc(q)+1.0&f>fc(q)-1.0)=nan;end
    X = fillmissing(X','nearest')';
  
    RatioA = sum(X(:,f<20 & f>1),2)./sum(X(:,f<55& f>1),2);
    RatioB = sum(X(:,f<9&f>7),2)./sum(X(:,f<10 & f>2),2);
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));

    %xState = ALL(:,idx);
    xState =  all_state;
    xState(xState==2)=1;
    Ss =xState;

    dS = Ss-[Ss(1);Ss(1:end-1)];
    beta = [RatioA(:)';RatioB(:)']; beta = fillmissing(beta,'nearest');
    
    T = find(dS==+2); % state 1
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(2,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==1&unique(Ss(T(tl)+(1:CutBins_2nd)))==3
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    ECoG_ALL_S1(:,:,idx)=nanmean(Gc,3);
    
    T = find(dS==-2); % state 2
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(2,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    ECoG_ALL_S2(:,:,idx)=nanmean(Gc,3);
    
    
    
end
x = squeeze(flip(sqrt(sum(BOLD_ALL_S1(1:4,:,:).^2,1)),2));
y = squeeze(    (sqrt(sum(BOLD_ALL_S2(1:4,:,:).^2,1))  ));
BOLD_Distance = abs((x-y)); BOLD_Distance=[BOLD_Distance(3:end,:);BOLD_Distance(1:2,:)];
x = squeeze(flip(sqrt(sum(ECoG_ALL_S1(1:2,:,:).^2,1)),2));
y = squeeze(    (sqrt(sum(ECoG_ALL_S2(1:2,:,:).^2,1))  ));
ECoG_Distance = abs((x-y)); ECoG_Distance=[ECoG_Distance(3:end,:);ECoG_Distance(1:2,:)];

F = figure; 
Ms = nanmean(BOLD_Distance,2);
SEMs = nanstd(BOLD_Distance,0,2)./sqrt(46);
t1=2*(-15:15);plot(t1,Ms); hold on;
patch('XData',[t1,fliplr(t1)],'YData',[Ms+SEMs;flipud(Ms-SEMs)]','facealpha',0.2,'edgecolor','none');
yyaxis right;
Ms = nanmean(ECoG_Distance,2);
SEMs = nanstd(ECoG_Distance,0,2)./sqrt(46);
t2=-30:30;plot(t2,Ms); hold on;
patch('XData',[t2,fliplr(t2)],'YData',[Ms+SEMs;flipud(Ms-SEMs)]','facealpha',0.2,'edgecolor','none');
ylim([0.0 0.04]);
xlim([-30 30])
set(gca,'tickdir','out')

x=BOLD_Distance;h1=zeros(size(x,1),1);for i=1:size(x,1);[h1(i)]=ttest(x(i,:),mean(x(1:10,:),1));end
x=ECoG_Distance;h2=zeros(size(x,1),1);for i=1:size(x,1);[h2(i)]=ttest(x(i,:),mean(x(1:20,:),1));end
plot(t1(h1==1),h1(h1==1)*0.02,'o');
plot(t2(h2==1),h2(h2==1)*0.03,'o');