clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\pca\';




load(fullfile('H:\ECoG\dcc_yyl\','states_tr.mat'));
ns = 300;xlc = 0;

%{
for shf = 00:1000
    
    
Sy=[];St=[];
for idx=1:46
    State = ALL(:,idx);
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
    State = cat(1,Ss{:});
    
    dS = State-[State(1);State(1:end-1)];
    x=load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta=x.beta;
    beta(:,numel(State)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(3,.1);
    beta = filtfilt(b,a,beta')';
    P = double(dS~=0);
    
    [~,LOCS] = findpeaks(P,'MinPeakHeight',0.5);
    LOCS=[LOCS;numel(P)+1];
    
    for xl=1:numel(LOCS)-1
        
        xlc=xlc+1;
        S = beta(1:10,LOCS(xl):LOCS(xl+1)-1);
        Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
        Si = fillmissing(Si','nearest')';
        Sy = cat(3,Sy,Si);
        
        
        S = State(LOCS(xl):LOCS(xl+1)-1)';
        Si = interp1((0:size(S,2)-1)/size(S,2)*ns, S', 1:ns, 'nearest')';
        Si = fillmissing(Si','nearest')';
        Si = median(S)*ones(size(Si));
        if median(S)==2
           1; 
        end
        St = cat(2,St,Si(:));
        
    end
end

dest = fullfile(WholePath,'StateStable');
if ~exist(dest,'dir');mkdir(dest);end

ty = median(St,1);
AVG1 = nanmean(Sy(:,:,ty==1),3);
AVG3 = nanmean(Sy(:,:,ty==3),3);
AVG4 = nanmean(Sy(:,:,ty==4),3);
save(fullfile(dest,['Shuffle_',num2str(shf,'%04d'),'.mat']),'AVG1','AVG3','AVG4');

if shf==0
    AVG1 = (Sy(:,:,ty==1));
    AVG3 = (Sy(:,:,ty==3));
    AVG4 = (Sy(:,:,ty==4));
    save(fullfile(dest,['Shuffle_',num2str(shf,'%04d'),'.mat']),'AVG1','AVG3','AVG4');
end

end
%}


clear Shuffle ;
for shf = 00:1000
    filepath = fullfile(WholePath,'StateStable');
    x = load(fullfile(filepath,['Shuffle_',num2str(shf,'%04d'),'.mat']));
    
    if shf==0; Raw = x; end
    if shf~=0; Shuffle(shf) = x; end
end


close all
for st = [1,3,4];
    
    F = figure;
    for pc=1:4
        eval(['x1=Raw.AVG',num2str(st),'(',num2str(pc),',:,:);']);
        eval(['y1=cat(1,Shuffle.AVG',num2str(st),');']);
        y2 = y1(pc:10:end,:);
        
        subplot(2,2,pc);
        x = squeeze(x1);
        Ms = mean(x,2); 
        SEMs = std(x,0,2)./sqrt(size(x,2));
        t = (1:300)';
        if pc==1;ColorVec=[217,83,25]/255;end
        if pc==2;ColorVec=[0 ,104,60]/255;end
        if pc==3;ColorVec=[140,182,44]/255;end
        if pc==4;ColorVec=[171,131,71]/255;end
        p=polarplot(t/max(t)*2*pi,Ms,'color',[ColorVec],'linewidth',2);
        hold on;
        theta = linspace(0,2*pi,300);
        RAM=[Ms(:)',mean(y2,1)-2.0*std(y2,0,1),mean(y2,1)+2.0*std(y2,0,1)];
        polarplot(t/max(t)*2*pi,mean(y2,1)-2.0*std(y2,0,1),'k-');
        polarplot(t/max(t)*2*pi,mean(y2,1)+2.0*std(y2,0,1),'k-');
        if min(RAM)<0;xlix(1)=min(RAM)*1.4;else;xlix(1)=min(RAM)*0.5;end
        if max(RAM)>0;xlix(2)=max(RAM)*1.4;else;xlix(2)=max(RAM)*0.5;end
        rlim([xlix(1) xlix(2)])
%         p=polarplot(theta,Ms,'color',[ColorVec],'linewidth',2);
        hold on
%         polarfill(gca,theta,mean(y2,1)-2.0*std(y2,0,1),mean(y2,1)+2.0*std(y2,0,1),[0 0 0],0.2);
%         hold on;
    end
    
    saveas(F,fullfile('H:\ECoG\pca',['State_',num2str(st),'.emf']));
    
    
    
    F = figure('Position', [592 762 1231 216]);
    for pc=1:4
        eval(['x1=Raw.AVG',num2str(st),'(',num2str(pc),',:,:);']);
        eval(['y1=cat(1,Shuffle.AVG',num2str(st),');']);
        y2 = y1(pc:10:end,:);
        
        subplot(1,4,pc);
        x = squeeze(x1);
        Ms = mean(x,2); 
        t = (1:300)';t0=t/max(t)*100;
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
    end
    saveas(F,fullfile('H:\ECoG\pca',['State_',num2str(st),'_control.emf']));
end


