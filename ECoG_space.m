clc;clear
close all
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\yyl_power\';

load(fullfile('H:\ECoG\dcc_yyl\','states_1s.mat'));

ALLState=[];
ALLRatioA=[];
ALLRatioB=[];

for idx = [1:15 17:46]
    
    
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
    X = fillmissing(X','linear')';

  
    RatioA = sum(X(:,f<20 & f>1),2)./sum(X(:,f<55& f>1),2);
    RatioB = sum(X(:,f<9&f>7),2)./sum(X(:,f<10 & f>2),2);
    
    ALLState = [ALLState;all_state(1:14400,1)];
    ALLRatioA = [ALLRatioA;RatioA(1:14400,1)];
    ALLRatioB = [ALLRatioB;RatioB(1:14400,1)];
    
    1;
    
end


State=ALLState;
State(State==2)=1;
uniqueS = unique(State);

dx = State-[State(1);State(1:end-1)];
P = double(dx~=0);
[~,LOCS] = findpeaks(P,'MinPeakHeight',0.5);
LOCS=[LOCS;numel(P)+1];
Sa=[];Sb=[];St=[]; ns = 100;
for xl=1:numel(LOCS)-1

    S = ALLRatioA(LOCS(xl):LOCS(xl+1)-1)';
    Si = interp1((0:size(S,2)-1)/size(S,2)*ns, S', 1:ns, 'nearest')';
    Sa(xl,:) = fillmissing(Si','nearest')';
    
    S = ALLRatioB(LOCS(xl):LOCS(xl+1)-1)';
    Si = interp1((0:size(S,2)-1)/size(S,2)*ns, S', 1:ns, 'nearest')';
    Sb(xl,:) = fillmissing(Si','nearest')';
    
    S = State(LOCS(xl):LOCS(xl+1)-1)';
    Si = interp1((0:size(S,2)-1)/size(S,2)*ns, S', 1:ns, 'nearest')';
    Si = fillmissing(Si','nearest')';
    St(xl,:) = median(S)*ones(size(Si));    
end



% close all
F = figure('Position',[680 511 497 467]);
for st = [1 3 4]
    RA = ALLRatioA(ALLState==st,:);
    RB = ALLRatioB(ALLState==st,:);
    NUM=30; 
    if st==4;NUM=30;end
    Rs=floor(size(RA,1)/NUM);
    QA = mean(reshape(RA(1:NUM*Rs,1),Rs,NUM),2);
    QB = mean(reshape(RB(1:NUM*Rs,1),Rs,NUM),2);
    if st==1;ColorVec=[232 176 097]/255;end
    if st==3;ColorVec=[016 129 054]/255;end
    if st==4;ColorVec=[017 037 103]/255;end
    
    if st~=4
        plot(QB(1:5:end),QA(1:5:end),'.','color',ColorVec);
    else
        plot(QB,QA,'.','color',ColorVec);
    end
    hold on;
    
    ty = median(St,2);
    y = nanmean(Sa(ty==st,:),1);
    x = nanmean(Sb(ty==st,:),1);
    [b_,a_]=butter(3,0.1);
    x=filtfilt(b_,a_,x);
    y=filtfilt(b_,a_,y);
    
%     plot(x,y,'k','linewidth',2);
end




% F = figure('Position',[680 511 497 467]);

Ss=ALLState; Ss(Ss==2)=1;
dS = Ss-[Ss(1);Ss(1:end-1)];

CutBins = 10;
for st = [ +1 -2 -3]
    
    T = find(dS==st); % state 1
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Ga = zeros(CutBins*2+1,numel(T));
    Gb = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        Ga(:,tl) = ALLRatioA(T(tl)+(-CutBins:CutBins));
        Gb(:,tl) = ALLRatioB(T(tl)+(-CutBins:CutBins));
    end
    [b_,a_]=butter(3,0.1);
    x=nanmean(Gb,2);%x=filtfilt(b_,a_,x);
    y=nanmean(Ga,2);%y=filtfilt(b_,a_,y);
%     plot(x);
    plot(x,y,'k','linewidth',2);
    hold on;
%     plot(y)
end


xlabel('EEG power B');
ylabel('EEG power A');
set(gca,'box','off','linewidth',2,'fontsize',15,'tickdir','out','ticklength',[0.04 1]);


