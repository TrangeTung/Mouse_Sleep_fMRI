clc;clear;
cd('H:\LFP_HIPP\1121');%cd('G:\ecog_yyl\hz');%
X = dir('m*');
g=41;%37,41,43,28 %20,10/17,1
for ix=1:numel(X)
    cd(fullfile(X(ix).folder,X(ix).name));
    load(['s',num2str(g),'.mat']);%s
    load('fsEcg.mat');
    load('fs.mat');
    ch5=cell2mat(struct2cell(load('ch5.mat')));%%%%%
    data = interp1((1:numel(ch5)) / 24414,ch5,1/1024:1/1024:numel(ch5) / 24414 , 'linear');   
    %% filter
    fs_raw=24414;
    fs_new=1024;
    data2(abs(data2)>800)=nan;
    data2 = fillmissing(data2,'linear');
    N=length(data2);
    n=0:N-1;
    fs=1024;
    x=1/fs_new:1/fs_new:length(data2)/fs_new;
    
    wn=[10*2 16*2]/fs;
    [k,l]=butter(2,wn);
    fhipp=filtfilt(k,l,data2);

    wn=[5*2 8*2]/fs;
    [k,l]=butter(2,wn);
    theta=filtfilt(k,l,data2);
%     
%     wn=[0.2*2 300*2]/fs;
%     [k,l]=butter(2,wn);
%     raw=filtfilt(k,l,sec_ch5);
    
    %% extract events
    env=envelope(fhipp,1000);
    env1=envelope(theta,1000);
    m=mean(env);
    threshold=mean(env)+1.5*std(env);
    threshold1=mean(env1)+1.5*std(env1);
    
    x01=(env>threshold);% spindle threshold
    x02=(env1>threshold);% theta threshold
    SE = strel('line',round(fs_new*0.8),0);
    x01_ = imdilate(x01,SE);
    x02_ = imdilate(x02,SE);
    xmerge = x01_(:).*x02_(:);
    SE = strel('line',round(fs_new*0.6),90); % exclude the short "spindle"
    xmerge = imerode(xmerge,SE);
    xmerge = imdilate(xmerge,SE);
    spindle1 = xmerge;
    
%     F = figure('Position', [30 558 1849 420]);
%     plot(x(1:5:end),sec_ch5(1:5:end)/10+80);
%     hold on;plot([min(x),max(x)],[110,110]);
%     plot(x(1:5:end),fhipp(1:5:end)/5);
%     plot(x(1:5:end),theta(1:5:end)/5+50);
%     plot(x(1:5:end),double(nrem_peak(1:5:end)*30)+25,'r','linewidth',2);
%     axis tight
%     ylim([-300 500]/2)
    
%% find all peaks & 去掉低于阈值的peak，剩余peak标为2
      
    [P,L]=findpeaks(env);
    peak=ones(1,14745600);
    peak(L)=2;
    for k=1:length(L)
        y=L(k);
        if env(y) < threshold
            peak(y)=1;
        else
            k=k+1;
        end
    end
    
    all_peak=spindle1.*peak';
    all_peak(14399*1024:14400*1024,1)=0;
%     figure;plot(all_peak*100);hold on;plot(env);plot(rip1);
   
%% 以peak为中心，寻找起始点和终点，并删除
    p=find(all_peak==2);
    for i=1:round(m)
        if p(i,1)<m
            p(i,1)=0;
        end
    end
    p(find(p==0))=[];

    
    spindle=zeros(1,14745600);
    n=1;
    for i=1:length(p)
        x=1;
        for x=1:p(i,1)
            bg=p(i,1)-x;
            if env(1,bg)<m || bg==1
                break
            else
                x=x+1;
            end
        end
        
        x=1;
        for x=1:p(i,1)
            ed=p(i,1)+x;
            if env(1,ed)<m
                break
            else
                x=x+1;
            end
        end
        
        if (ed-bg)<0.4*1024 || (ed-bg)>3*1024
        else
            piece(n,1)=bg;
            piece(n,3)=ed;
            spindle(1,bg:ed)=1;
            n=n+1;
        end
    end
    
%     figure;plot(ripple*30+30);hold on;plot(rip1);plot(env);
    
    all_cut=unique(piece,'rows');
    
    for i=1:length(all_cut)
        f1=all_cut(i,1);
        f2=all_cut(i,3);
        [Y,I]=max(fhipp(1,f1:f2));
        all_cut(i,2)=I+f1;
    end
    
    spindle(all_cut(:,2))=2;
%     figure;plot(ripple*30+30);hold on;plot(rip1/5);plot(env/5);
%% brainstate

    for i=1:14400
        all_state1(((i-1)*1024+1):((i-1)*1024+1024),1)=all_state(i,1);
    end
   
    e=zeros(14745600,1);
    e1=zeros(14745600,1);
    e(find(all_state1==1))=1;
    e1(find(all_state1==3))=1;
    awake_peak=spindle.*e';
    nrem_peak=spindle.*e1';
    
%     n=zeros(1,14745600);
%     n(1,:)=threshold;
%     q=zeros(1,14745600);
%     q(1,:)=threshold1;
%     figure;plot(sec_ch5(1,7578300:7588540)/5+200);hold on;plot(env(1,7578300:7588540)/2,'k');plot(fhipp(1,7578300:7588540)/2,'k');plot(n(1,7578300:7588540)/2,'--');plot(nrem_peak(1,7578300:7588540)*30+500);axis tight;%plot(theta(1,7579000:7589240)/2-200,'k');plot(q(1,7579000:7589240)/2-200,'--');plot(env1(1,7579000:7589240)/2-200,'k');
%     figure;plot(sec_ch5(1,7583000:7586072)/5+300);hold on;plot(env(1,7583000:7586072)/2,'k');plot(fhipp(1,7583000:7586072)/2,'k');plot(n(1,7583000:7586072)/2,'--');plot(nrem_peak(1,7583000:7586072)*30+500);axis tight;

    
    awake_cut=[];
    nrem_cut=[];
    for i=1:length(all_cut)
        if awake_peak(all_cut(i,1))==1
            awake_cut=[awake_cut;all_cut(i,:)];
        else
            nrem_cut=[nrem_cut;all_cut(i,:)];
        end
    end

    figure;plot(fhipp);hold on;plot(nrem_peak*30+50,'r','linewidth',2);axis tight
    spindle = nrem_peak;
    cut = nrem_cut;
    
    save(['spindle_peak',num2str(g),'.mat'],'spindle','cut');
    g=g+1;
end

